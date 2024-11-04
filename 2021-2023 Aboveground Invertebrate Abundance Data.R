#### Libraries ####
library(readxl)
library(tidyverse)
library(codyn)
library(ggplot2)  
library(ggpubr)
library(gridExtra)
library(vegan)
library(nlme)
library(cowplot)
library(emmeans)
library(openxlsx)
#### Seed Set ####
set.seed(123)

#### CSV read ####

#Contains aboveground Abundance ID information
Abundance_ID_aboveground <- read_excel("Aboveground Invertebrate ID.xlsx")
Abundance_ID_aboveground <- as_tibble(Abundance_ID_aboveground)

#Contains time sense burned information for each watershed
BurnInfo2021 <- read_excel("YearsSenseBurned2021.xlsx")
BurnInfo2022 <- read_excel("YearsSenseBurned2022.xlsx")
BurnInfo2023 <- read_excel("YearsSenseBurned2023.xlsx")

Weight_2 <- read_excel("Aboveground Invertebrate Biomass.xlsx")


#### Aboveground Abundance Data Clean Up ####

#Checking for WS typo
table(Abundance_ID_aboveground$WS)

#Fixing typos
Abundance_ID_aboveground$WS <- gsub("C35A", "C3SA", Abundance_ID_aboveground$WS)
Abundance_ID_aboveground$WS <- gsub("C3D", "C3B", Abundance_ID_aboveground$WS)
Abundance_ID_aboveground$WS <- gsub("CIA", "C1A", Abundance_ID_aboveground$WS)

#Typo fix confirmation
table(Abundance_ID_aboveground$WS)


#### Above ground Weight Data Clean up ####

#Subtracting tube weight from final mass to get sample weight. Also setting negatives to zero (they are true zeros, the scale just doesn't register)
Weight_2$weight <- as.numeric(Weight_2$`Final mass (mg)`) - as.numeric(Weight_2$`Tube mass (mg)`)
Weight_2$weight[Weight_2$weight < 0] <- 0

#Create a sample column with Date, Watershed, Transect and Plot
Weight_2$sample <- paste(Weight_2$Date, Weight_2$Watershed, Weight_2$Transect, Weight_2$Plot, sep = "_")

Weight_2 <- Weight_2[!grepl("HUGE GRASSSHOPPER V", Weight_2$Notes),]


Weight_2_combined <- Weight_2 %>%
  group_by(sample) %>%
  summarise(combined_weight = sum(weight, na.rm = TRUE))

#Fix typo
Weight_2_combined$sample <- gsub("2012_C3A_A_1", "2021_C3A_A_1", Weight_2_combined$sample)

# Adding useful information
Weight_2_combined <- Weight_2_combined %>%
  mutate(Treatment = ifelse(grepl("C1", sample), "ABG", "PBG"),
         block = ifelse(grepl("S", sample), "North", "South"))



#### Getting Community Metrics ####


#Getting  community data and tagging ABG and PBG

#Make a new Count column. Solves weird non-numeric column problems.

# Convert "Observed #" column to numeric
Abundance_ID_aboveground$`Observed #` <- as.numeric(Abundance_ID_aboveground$`Observed #`)

# Calculate the sum of "Observed #" and "Collected #" columns
Abundance_ID_aboveground$Count <- rowSums(select(Abundance_ID_aboveground, `Observed #`, `Collected #`), na.rm = TRUE)

# Replace 0s or NAs with 1s in the "Count" column
Abundance_ID_aboveground$Count <- ifelse(Abundance_ID_aboveground$Count == 0 | is.na(Abundance_ID_aboveground$Count), 1, Abundance_ID_aboveground$Count)

Abundance_ID_aboveground <- Abundance_ID_aboveground %>%
  mutate(
    # Convert "Observed #" column to numeric
    `Observed #` = as.numeric(`Observed #`),
    # Calculate the sum of "Observed #" and "Collected #" columns
    Count = rowSums(select(., `Observed #`, `Collected #`), na.rm = TRUE),
    # Replace 0s or NAs with 1s in the "Count" column
    Count = ifelse(Count == 0 | is.na(Count), 1, Count)
  )

Abundance_ID_aboveground <- Abundance_ID_aboveground %>% 
  mutate("ID" = paste(Order, Family, sep = "_")) %>% 
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

# Joined full_join(BurnInfo, commMetrics2) %>% unite("TreatmentSB",c("Treatment","SB"), sep="_") 

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
  separate(Sample, into = c("WS", "Trans", "Plot"), sep = "_", extra = "drop") %>% mutate(block = ifelse(grepl("S", WS), "North", "South"))



counts <- ggplot(total_counts, aes(x = Count, color = interaction(Treatment, block), linetype = interaction(Treatment, block))) +
  geom_density() +
  labs(title = "Aboveground Invertebrate Count Density Plot",
       x = "Count",
       y = "Density") +
  scale_color_manual(values = rep(c("blue", "red", "blue", "red"), 2)) +
  scale_linetype_manual(values = rep(c("solid", "dashed"), each = 2)) +
  theme_classic()

print(counts)

Weight <- ggplot(Weight_2_combined, aes(x = combined_weight, color = interaction(Treatment, block), linetype = interaction(Treatment, block))) +
  geom_density() +
  labs(title = "Aboveground Invertebrate Weight Density Plot",
       x = "Weight",
       y = "Density") +
  scale_color_manual(values = rep(c("blue", "red", "blue", "red"), 2)) +
  scale_linetype_manual(values = rep(c("solid", "dashed"), each = 2)) +
  theme_classic()

print(Weight)

grid.arrange(richness, evenness, counts, Weight)

#Without weight

grid.arrange(richness, evenness, counts, nrow =2, ncol = 2)




#### Graph Setting for NMDS  ####
NMDS_AXIS_TEXT_SIZE <- 50
STATS_TEXT_SIZE <- 18
#### 2021 PERMANOVA & NMDS ####
### by Years since burn

#Combing count with burn info


#merge with burn data
Abundance_Data2021 <- full_join(Abundance_ID_aboveground, BurnInfo2021, by = "WS") %>% 
  unite("TreatmentSB",c("Treatment","SB"), sep="_") %>% 
  filter(Date == "2021") %>%   
  filter(Plot %in% c("2", "4"))




#tag abg and pbg
Abundance_Data2021$Sample <- paste(Abundance_Data2021$WS, Abundance_Data2021$Trans, Abundance_Data2021$Plot, sep = "_")
Abundance_Data2021$Treatment <- ifelse(grepl("C1", Abundance_Data2021$Sample), "ABG", "PBG") 

# PERMANOVA


abundanceWide <- Abundance_Data2021 %>% 
  mutate(block = ifelse(grepl("S", Sample), "North", "South")) %>%
  select(Sample, block, WS, Trans, Plot, Treatment, TreatmentSB, Count, ID) %>%
  group_by(Sample, block, WS, Trans, Plot, Treatment, TreatmentSB, ID) %>% 
  summarise(Count = sum(Count)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = 'ID', values_from = 'Count', values_fill = list(Count = 0)) 

abundanceWide <- abundanceWide %>% 
  mutate(sum = rowSums(abundanceWide[, c(8:150)], na.rm = TRUE)) %>% 
  filter(sum > 0, !(Sample %in% c('C1A_A_38_ABG', 'C3SA_D_16_PBG', 'C3SA_C_38_PBG')))

abundanceWide_2021 <- abundanceWide

print(permanova <- adonis2(formula = abundanceWide[,8:150]~TreatmentSB, data=abundanceWide, permutations=999, method="bray"))
#F=0.7814, df=3,57, p=0.86

#pairwise

pairwise_results <- pairwise.adonis(abundanceWide[, 8:150], abundanceWide$TreatmentSB)


#betadisper
veg <- vegdist(abundanceWide[,8:150], method = "bray")
dispersion <- betadisper(veg, abundanceWide$TreatmentSB)
permutest(dispersion, pairwise=TRUE, permutations=999) 
#F=0.5209, df=3,57, p=0.663 

BC_Data <- metaMDS(abundanceWide[,8:150])
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
Years_Since_Burned_NMDS_2021 <- ggplot(BC_NMDS_Graph, aes(x = MDS1, y = MDS2, color = group, linetype = group, shape = group)) +
  geom_point(size = 6) + 
  geom_path(data = BC_Ellipses, aes(x = NMDS1, y = NMDS2), size = 3) +
  labs(color = "", linetype = "", shape = "") +
  scale_colour_manual(values = c("blue", "#7A2021", "#C21A09", "#FF0800"), name = "",
                      labels = c('ABG', 'PBG 0', 'PBG 1', 'PBG 2')) +
  scale_linetype_manual(values = c("solid", "twodash", "twodash", "twodash"), name = "",
                        labels = c('ABG', 'PBG 0', 'PBG 1', 'PBG 2')) +
  scale_shape_manual(values = c(19, 19, 17, 15),
                     labels = c('ABG', 'PBG 0', 'PBG 1', 'PBG 2')) +
  xlab("NMDS1") + 
  ylab("NMDS2") + 
  theme_bw() +
  theme(axis.text.x = element_text(size = NMDS_AXIS_TEXT_SIZE, color = "black"), 
        axis.text.y = element_text(size = NMDS_AXIS_TEXT_SIZE, color = "black"), 
        axis.title.x = element_text(size = NMDS_AXIS_TEXT_SIZE, color = 'black'),
        axis.title.y = element_text(size = NMDS_AXIS_TEXT_SIZE, color = 'black'),
        legend.text = element_text(size = 60),
        legend.key.size = unit(5, "cm"),  # Adjust the size of legend key (dots)
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        plot.title = element_text(size = 75)
  ) +
  ggtitle("2021 year since burned") + 
  annotate('text', x = min(BC_NMDS_Graph$MDS1) - 0.04, y = min(BC_NMDS_Graph$MDS2) + 0.06, 
           label = 'Mean p = 0.069\nVariance p = 0.451', size = STATS_TEXT_SIZE, hjust = 'left')

Years_Since_Burned_NMDS_2021

Years_Since_Burned_NMDS_2021

#PERMANOVA: F=1.2983  , df=3,57, p=0.069 
#Betadisperion: #F=0.8876    , df=3,57, p=0.451 


#export at 1500x1000

### by watershed

Abundance_Data2021 <- full_join(Abundance_ID_aboveground, BurnInfo2021, by = "WS") %>% 
  unite("TreatmentSB",c("Treatment","SB"), sep="_") %>% 
  filter(Date == "2021") %>%   
  filter(Plot %in% c("2", "4"))

#Filter ABG 
ABG_Test <- full_join(Abundance_ID_aboveground, BurnInfo2021, by = "WS") %>% 
  unite("TreatmentSB",c("Treatment","SB"), sep="_") %>% 
  filter(Date == "2021") %>%   
  filter(Plot %in% c("2", "4")) %>% 
  filter(TreatmentSB == "ABG_0")

#Filter PBG
PBG_Test <- full_join(Abundance_ID_aboveground, BurnInfo2021, by = "WS") %>% 
  unite("TreatmentSB",c("Treatment","SB"), sep="_") %>% 
  filter(Date == "2021") %>%   
  filter(Plot %in% c("2", "4")) %>% 
  filter(TreatmentSB == c("PBG_0", "PBG_1", "PBG_2"))


# Set seed for reproducibility
set.seed(123)

# Get unique samples
unique_samples <- unique(PBG_Test$Sample)

# Randomly select 16 unique samples
subsamples <- sample(unique_samples, 16, replace = FALSE)

# Filter the data frame based on the selected samples
subsampled_data <- PBG_Test %>% filter(Sample %in% subsamples)

#New Abundance_Data2021 with subsamples
Abundance_Data2021 <- full_join(subsampled_data, ABG_Test)

#tag abg and pbg
Abundance_Data2021$Sample <- paste(Abundance_Data2021$WS, Abundance_Data2021$Trans, Abundance_Data2021$Plot, sep = "_")
Abundance_Data2021$Treatment <- ifelse(grepl("C1", Abundance_Data2021$Sample), "ABG", "PBG") 

# PERMANOVA


abundanceWide <- Abundance_Data2021 %>% 
  mutate(block = ifelse(grepl("S", Sample), "North", "South")) %>%
  select(Sample, block, WS, Trans, Plot, Treatment, TreatmentSB, Count, ID) %>%
  group_by(Sample, block, WS, Trans, Plot, Treatment, TreatmentSB, ID) %>% 
  summarise(Count = sum(Count)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = 'ID', values_from = 'Count', values_fill = list(Count = 0)) 

abundanceWide <- abundanceWide %>% 
  mutate(sum = rowSums(abundanceWide[, c(8:89)], na.rm = TRUE)) %>% 
  filter(sum > 0, !(Sample %in% c('C1A_A_38_ABG', 'C3SA_D_16_PBG', 'C3SA_C_38_PBG')))


# PERMANOVA
print(permanova <- adonis2(formula = abundanceWide[,8:89]~Treatment, data=abundanceWide, permutations=999, method="bray"))
#F=1.9529 , df=1,28, p=0.013 *

#betadisper
veg <- vegdist(abundanceWide[,8:89], method = "bray")
dispersion <- betadisper(veg, abundanceWide$Treatment)
permutest(dispersion, pairwise=TRUE, permutations=999) 
#F=23.434       , df=1,28, p=0.001 ***

BC_Data <- metaMDS(abundanceWide[,8:89])
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
ABG_VS_PBG_NMDS_2021 <- ggplot(BC_NMDS_Graph, aes(x = MDS1, y = MDS2, color = group, linetype = group, shape = group)) +
  geom_point(size = 6) + 
  geom_path(data = BC_Ellipses, aes(x = NMDS1, y = NMDS2), size = 3) +
  labs(color = "", linetype = "", shape = "") +
  scale_colour_manual(values = c("blue", "red"), name = "") +
  scale_linetype_manual(values = c("solid", "twodash"), name = "") +
  scale_shape_manual(values = c(19, 19)) +
  xlab("NMDS1") + 
  ylab("NMDS2") + 
  annotate('text', x = min(BC_NMDS_Graph$MDS1) - 0.04, y = min(BC_NMDS_Graph$MDS2) + 0.06, 
           label = 'Mean p = 0.013\nVariance p = 0.001', size = STATS_TEXT_SIZE, hjust = 'left') +
  theme_classic() +
  theme(axis.text.x = element_text(size = NMDS_AXIS_TEXT_SIZE, color = "black"), 
        axis.text.y = element_text(size = NMDS_AXIS_TEXT_SIZE, color = "black"), 
        axis.title.x = element_text(size = NMDS_AXIS_TEXT_SIZE, color = 'black'),
        axis.title.y = element_text(size = NMDS_AXIS_TEXT_SIZE, color = 'black'),
        legend.text = element_text(size = 60),
        legend.key.size = unit(5, "cm"),  # Adjust the size of legend key (dots)
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        plot.title = element_text(size = 75)
  ) +
  ggtitle("2021 ABG vs PBG")


#PERMANOVA: F=1.9529 , df=1,28, p=0.013 *
#PERMUTEST: F=23.434       , df=1,28, p=0.001 ***

ABG_VS_PBG_NMDS_2021


#### 2022 PERMANOVA & NMDS ####
# Combine count with burn info


# Merge with burn data
Abundance_Data2022 <- full_join(Abundance_ID_aboveground, BurnInfo2022, by = "WS") %>% 
  unite("TreatmentSB", c("Treatment", "SB"), sep = "_") %>% 
  filter(Date == "2022")


# Tag ABG and PBG
Abundance_Data2022$Sample <- paste(Abundance_Data2022$WS, Abundance_Data2022$Trans, Abundance_Data2022$Plot, sep = "_")
Abundance_Data2022$Treatment <- ifelse(grepl("C1", Abundance_Data2022$Sample), "ABG", "PBG") 

# PERMANOVA
abundanceWide <- Abundance_Data2022 %>% 
  mutate(block = ifelse(grepl("S", Sample), "North", "South")) %>%
  select(Sample, block, WS, Trans, Plot, Treatment, TreatmentSB, Count, ID) %>%
  group_by(Sample, block, WS, Trans, Plot, Treatment, TreatmentSB, ID) %>% 
  summarise(Count = sum(Count)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = 'ID', values_from = 'Count', values_fill = list(Count = 0)) 

abundanceWide <- abundanceWide %>% 
  mutate(sum = rowSums(abundanceWide[, c(8:79)], na.rm = TRUE)) %>% 
  filter(sum > 0, !(Sample %in% c('C1A_A_38_ABG', 'C3SA_D_16_PBG', 'C3SA_C_38_PBG')))

abundanceWide_2022 <- abundanceWide


print(permanova <- adonis2(formula = abundanceWide[, 8:121] ~ TreatmentSB, data = abundanceWide, permutations = 999, method = "bray"))
# F=2.1065, df=3,60, p=0.001 ***



# Betadisper
veg <- vegdist(abundanceWide[, 8:121], method = "bray")
dispersion <- betadisper(veg, abundanceWide$TreatmentSB)
permutest(dispersion, pairwise = TRUE, permutations = 999) 
# F=0.0305, df=3,60, p=0.992

BC_Data <- metaMDS(abundanceWide[, 8:121])
sites <- 1:nrow(abundanceWide)
BC_Meta_Data <- abundanceWide[, 1:7]
plot(BC_Data$points, col = as.factor(BC_Meta_Data$TreatmentSB))
ordiellipse(BC_Data, groups = as.factor(BC_Meta_Data$TreatmentSB), kind = "sd", display = "sites", label = TRUE)

veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) {
  theta <- (0:npoints) * 2 * pi / npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

# Generate ellipses
BC_NMDS <- data.frame(MDS1 = BC_Data$points[, 1], MDS2 = BC_Data$points[, 2], group = BC_Meta_Data$TreatmentSB)
BC_NMDS_Graph <- cbind(BC_Meta_Data, BC_NMDS)
BC_Ord_Ellipses <- ordiellipse(BC_Data, BC_Meta_Data$TreatmentSB, display = "sites", kind = "se", conf = 0.95, label = TRUE)

BC_Ellipses <- data.frame()
# Generate ellipses points
for (g in levels(as.factor(BC_NMDS$group))) {
  BC_Ellipses <- rbind(BC_Ellipses,
                       cbind(as.data.frame(with(BC_NMDS[BC_NMDS$group == g, ], 
                                                veganCovEllipse(BC_Ord_Ellipses[[g]]$cov,
                                                                BC_Ord_Ellipses[[g]]$center,
                                                                BC_Ord_Ellipses[[g]]$scale))),
                             group = g))
}

Years_Since_Burned_NMDS_2022 <- ggplot(BC_NMDS_Graph, aes(x = MDS1, y = MDS2, color = group, linetype = group, shape = group)) +
  geom_point(size = 6) + 
  geom_path(data = BC_Ellipses, aes(x = NMDS1, y = NMDS2), size = 3) +
  labs(color = "", linetype = "", shape = "") +
  scale_colour_manual(values = c("blue", "#7A2021", "#C21A09", "#FF0800"), name = "",
                      labels = c('ABG', 'PBG 0', 'PBG 1', 'PBG 2')) +
  scale_linetype_manual(values = c("solid", "twodash", "twodash", "twodash"), name = "",
                        labels = c('ABG', 'PBG 0', 'PBG 1', 'PBG 2')) +
  scale_shape_manual(values = c(19, 19, 17, 15),
                     labels = c('ABG', 'PBG 0', 'PBG 1', 'PBG 2')) +
  xlab("NMDS1") + 
  ylab("NMDS2") + 
  theme_bw() +
  theme(axis.text.x = element_text(size = NMDS_AXIS_TEXT_SIZE, color = "black"), 
        axis.text.y = element_text(size = NMDS_AXIS_TEXT_SIZE, color = "black"), 
        axis.title.x = element_text(size = NMDS_AXIS_TEXT_SIZE, color = 'black'),
        axis.title.y = element_text(size = NMDS_AXIS_TEXT_SIZE, color = 'black'),
        legend.text = element_text(size = 60),
        legend.key.size = unit(5, "cm"),  # Adjust the size of legend key (dots)
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        plot.title = element_text(size = 75)
  ) +
  ggtitle("2022 year since burned") + 
  annotate('text', x = min(BC_NMDS_Graph$MDS1) - 0.04, y = min(BC_NMDS_Graph$MDS2) + 0.12, 
           label = 'Mean p = 0.001\nVariance p = 0.975', size = STATS_TEXT_SIZE, hjust = 'left')

# Mean F=2.4861, df=3,60, p=0.002 **
# Variance F=0.0305, df=3,60, p=0.992

Years_Since_Burned_NMDS_2022

### by watershed
#Filter ABG 
ABG_Test <- full_join(Abundance_ID_aboveground, BurnInfo2022, by = "WS") %>% 
  unite("TreatmentSB",c("Treatment","SB"), sep="_") %>% 
  filter(Date == "2022") %>%   
  filter(TreatmentSB == "ABG_0")

#Filter PBG
PBG_Test <- full_join(Abundance_ID_aboveground, BurnInfo2022, by = "WS") %>% 
  unite("TreatmentSB", c("Treatment", "SB"), sep = "_") %>% 
  filter(Date == "2022") %>%   
  filter(TreatmentSB %in% c("PBG_0", "PBG_1", "PBG_2"))



# Set seed for reproducibility
set.seed(123)

# Get unique samples
unique_samples <- unique(PBG_Test$Sample)

# Randomly select 16 unique samples
subsamples <- sample(unique_samples, 16, replace = FALSE)

# Filter the data frame based on the selected samples
subsampled_data <- PBG_Test %>% filter(Sample %in% subsamples)

#New Abundance_Data2021 with subsamples
Abundance_Data2022 <- full_join(subsampled_data, ABG_Test)


# Tag ABG and PBG
Abundance_Data2022$Sample <- paste(Abundance_Data2022$WS, Abundance_Data2022$Trans, Abundance_Data2022$Plot, sep = "_")
Abundance_Data2022$Treatment <- ifelse(grepl("C1", Abundance_Data2022$Sample), "ABG", "PBG") 

abundanceWide <- Abundance_Data2022 %>% 
  mutate(block = ifelse(grepl("S", Sample), "North", "South")) %>%
  select(Sample, block, WS, Trans, Plot, Treatment, TreatmentSB, Count, ID) %>%
  group_by(Sample, block, WS, Trans, Plot, Treatment, TreatmentSB, ID) %>% 
  summarise(Count = sum(Count)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = 'ID', values_from = 'Count', values_fill = list(Count = 0)) 

abundanceWide <- abundanceWide %>% 
  mutate(sum = rowSums(abundanceWide[, c(8:79)], na.rm = TRUE)) %>% 
  filter(sum > 0, !(Sample %in% c('C1A_A_38_ABG', 'C3SA_D_16_PBG', 'C3SA_C_38_PBG')))


# PERMANOVA
abundanceWide <- Abundance_Data2022 %>% 
  mutate(block = ifelse(grepl("S", Sample), "North", "South")) %>%
  select(Sample, block, WS, Trans, Plot, Treatment, TreatmentSB, Count, ID) %>%
  group_by(Sample, block, WS, Trans, Plot, Treatment, TreatmentSB, ID) %>% 
  summarise(Count = sum(Count)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = 'ID', values_from = 'Count', values_fill = list(Count = 0)) 

abundanceWide <- abundanceWide %>% 
  mutate(sum = rowSums(abundanceWide[, c(8:79)], na.rm = TRUE)) %>% 
  filter(sum > 0, !(Sample %in% c('C1A_A_38_ABG', 'C3SA_D_16_PBG', 'C3SA_C_38_PBG')))


print(permanova <- adonis2(formula = abundanceWide[,8:79]~Treatment, data=abundanceWide, permutations=999, method="bray"))
#F=4.5219, df=1,30, p=0.001 ***

#betadisper
veg <- vegdist(abundanceWide[,8:79], method = "bray")
dispersion <- betadisper(veg, abundanceWide$Treatment)
permutest(dispersion, pairwise=TRUE, permutations=999) 

#F=59.943 , df=1,30, p=0.001 ***

BC_Data <- metaMDS(abundanceWide[,8:79])
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
ABG_VS_PBG_NMDS_2022  <- ggplot(BC_NMDS_Graph, aes(x = MDS1, y = MDS2, color = group, linetype = group, shape = group)) +
  geom_point(size = 6) + 
  geom_path(data = BC_Ellipses, aes(x = NMDS1, y = NMDS2), size = 3) +
  labs(color = "", linetype = "", shape = "") +
  scale_colour_manual(values = c("blue", "red"), name = "") +
  scale_linetype_manual(values = c("solid", "twodash"), name = "") +
  scale_shape_manual(values = c(19, 19)) +
  xlab("NMDS1") + 
  ylab("NMDS2") +
  theme_classic() +
  theme(axis.text.x = element_text(size = NMDS_AXIS_TEXT_SIZE, color = "black"), 
        axis.text.y = element_text(size = NMDS_AXIS_TEXT_SIZE, color = "black"), 
        axis.title.x = element_text(size = NMDS_AXIS_TEXT_SIZE, color = 'black'),
        axis.title.y = element_text(size = NMDS_AXIS_TEXT_SIZE, color = 'black'),
        legend.text = element_text(size = 60),
        legend.key.size = unit(5, "cm"),  # Adjust the size of legend key (dots)
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        plot.title = element_text(size = 75)
  ) +
  ggtitle("2022 ABG vs PBG") +
  annotate('text', x = min(BC_NMDS_Graph$MDS1) - 0.04, y = min(BC_NMDS_Graph$MDS2) + 0.12, label = 'Variance p = 0.001', size = STATS_TEXT_SIZE, hjust = 'left') +
  annotate('text', x = min(BC_NMDS_Graph$MDS1) - 0.04, y = min(BC_NMDS_Graph$MDS2) + 0.010, label = 'Mean p=0.001', size = STATS_TEXT_SIZE, hjust = 'left')


#Mean F=4.5219, df=1,30, p=0.001 ***
#Variance F=59.943 , df=1,30, p=0.001 ***

ABG_VS_PBG_NMDS_2022


#### 2023 PERMANOVA & NMDS ####
# Combine count with burn info

# Merge with burn data
Abundance_Data2023 <- full_join(Abundance_ID_aboveground, BurnInfo2023, by = "WS") %>% 
  unite("TreatmentSB", c("Treatment", "SB"), sep = "_") %>% 
  filter(Date == "2023")

# Tag ABG and PBG
Abundance_Data2023$Sample <- paste(Abundance_Data2023$WS, Abundance_Data2023$Trans, Abundance_Data2023$Plot, sep = "_")
Abundance_Data2023$Treatment <- ifelse(grepl("C1", Abundance_Data2023$Sample), "ABG", "PBG") 

# PERMANOVA
abundanceWide <- Abundance_Data2023 %>% 
  mutate(block = ifelse(grepl("S", Sample), "North", "South")) %>%
  select(Sample, block, WS, Trans, Plot, Treatment, TreatmentSB, Count, ID) %>%
  group_by(Sample, block, WS, Trans, Plot, Treatment, TreatmentSB, ID) %>% 
  summarise(Count = sum(Count)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = 'ID', values_from = 'Count', values_fill = list(Count = 0)) 

abundanceWide <- abundanceWide %>% 
  mutate(sum = rowSums(abundanceWide[, c(8:75)], na.rm = TRUE)) %>% 
  filter(sum > 0, !(Sample %in% c('C1A_A_38_ABG', 'C3SA_D_16_PBG', 'C3SA_C_38_PBG')))

abundanceWide_2023  <- abundanceWide


print(permanova <- adonis2(formula = abundanceWide[, 8:75] ~ TreatmentSB, data = abundanceWide, permutations = 999, method = "bray"))
# F=1.3803, df=3,60, p=0.094 

# Betadisper
veg <- vegdist(abundanceWide[, 8:75], method = "bray")
dispersion <- betadisper(veg, abundanceWide$TreatmentSB)
permutest(dispersion, pairwise = TRUE, permutations = 999) 
# F=0.0039, df=3,60, p=0.999

BC_Data <- metaMDS(abundanceWide[, 8:75])
sites <- 1:nrow(abundanceWide)
BC_Meta_Data <- abundanceWide[, 1:7]
plot(BC_Data$points, col = as.factor(BC_Meta_Data$TreatmentSB))
ordiellipse(BC_Data, groups = as.factor(BC_Meta_Data$TreatmentSB), kind = "sd", display = "sites", label = TRUE)

veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) {
  theta <- (0:npoints) * 2 * pi / npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

# Generate ellipses
BC_NMDS <- data.frame(MDS1 = BC_Data$points[, 1], MDS2 = BC_Data$points[, 2], group = BC_Meta_Data$TreatmentSB)
BC_NMDS_Graph <- cbind(BC_Meta_Data, BC_NMDS)
BC_Ord_Ellipses <- ordiellipse(BC_Data, BC_Meta_Data$TreatmentSB, display = "sites", kind = "se", conf = 0.95, label = TRUE)

BC_Ellipses <- data.frame()
# Generate ellipses points
for (g in levels(as.factor(BC_NMDS$group))) {
  BC_Ellipses <- rbind(BC_Ellipses,
                       cbind(as.data.frame(with(BC_NMDS[BC_NMDS$group == g, ], 
                                                veganCovEllipse(BC_Ord_Ellipses[[g]]$cov,
                                                                BC_Ord_Ellipses[[g]]$center,
                                                                BC_Ord_Ellipses[[g]]$scale))),
                             group = g))
}

Years_Since_Burned_NMDS_2023 <- ggplot(BC_NMDS_Graph, aes(x = MDS1, y = MDS2, color = group, linetype = group, shape = group)) +
  geom_point(size = 6) + 
  geom_path(data = BC_Ellipses, aes(x = NMDS1, y = NMDS2), size = 3) +
  labs(color = "", linetype = "", shape = "") +
  scale_colour_manual(values = c("blue", "#7A2021", "#C21A09", "#FF0800"), name = "",
                      labels = c('ABG', 'PBG 0', 'PBG 1', 'PBG 2')) +
  scale_linetype_manual(values = c("solid", "twodash", "twodash", "twodash"), name = "",
                        labels = c('ABG', 'PBG 0', 'PBG 1', 'PBG 2')) +
  scale_shape_manual(values = c(19, 19, 17, 15),
                     labels = c('ABG', 'PBG 0', 'PBG 1', 'PBG 2')) +
  xlab("NMDS1") + 
  ylab("NMDS2") + 
  theme_bw() +
  theme(axis.text.x = element_text(size = NMDS_AXIS_TEXT_SIZE, color = "black"), 
        axis.text.y = element_text(size = NMDS_AXIS_TEXT_SIZE, color = "black"), 
        axis.title.x = element_text(size = NMDS_AXIS_TEXT_SIZE, color = 'black'),
        axis.title.y = element_text(size = NMDS_AXIS_TEXT_SIZE, color = 'black'),
        legend.text = element_text(size = 60),
        legend.key.size = unit(5, "cm"),  # Adjust the size of legend key (dots)
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        plot.title = element_text(size = 75)
  ) +
  ggtitle("2023 year since burned") + 
  annotate('text', x = min(BC_NMDS_Graph$MDS1) - 0.04, y = c(min(BC_NMDS_Graph$MDS2) + 0.03, min(BC_NMDS_Graph$MDS2) - 0.2), 
           label = c('Mean p = 0.094', 'Variance p = 0.999'), size = STATS_TEXT_SIZE, hjust = 'left')

# mean F=1.3803, df=3,60, p=0.094 
# variance F=0.0039, df=3,60, p=0.999

Years_Since_Burned_NMDS_2023

### by watershed

#Filter ABG 
ABG_Test <- full_join(Abundance_ID_aboveground, BurnInfo2023, by = "WS") %>% 
  unite("TreatmentSB",c("Treatment","SB"), sep="_") %>% 
  filter(Date == "2023") %>%   
  filter(TreatmentSB == "ABG_0")

#Filter PBG
PBG_Test <- full_join(Abundance_ID_aboveground, BurnInfo2023, by = "WS") %>% 
  unite("TreatmentSB",c("Treatment","SB"), sep="_") %>% 
  filter(Date == "2023") %>%   
  filter(TreatmentSB == c("PBG_0", "PBG_1", "PBG_2"))


# Set seed for reproducibility
set.seed(234)

# Get unique samples
unique_samples <- unique(PBG_Test$Sample)

# Randomly select 16 unique samples
subsamples <- sample(unique_samples, 16, replace = FALSE)

# Filter the data frame based on the selected samples
subsampled_data <- PBG_Test %>% filter(Sample %in% subsamples)

#New Abundance_Data2021 with subsamples
Abundance_Data2023 <- full_join(subsampled_data, ABG_Test)

# Tag ABG and PBG
Abundance_Data2023$Sample <- paste(Abundance_Data2023$WS, Abundance_Data2023$Trans, Abundance_Data2023$Plot, sep = "_")
Abundance_Data2023$Treatment <- ifelse(grepl("C1", Abundance_Data2023$Sample), "ABG", "PBG") 

# PERMANOVA
abundanceWide <- Abundance_Data2023 %>% 
  mutate(block = ifelse(grepl("S", Sample), "North", "South")) %>%
  select(Sample, block, WS, Trans, Plot, Treatment, TreatmentSB, Count, ID) %>%
  group_by(Sample, block, WS, Trans, Plot, Treatment, TreatmentSB, ID) %>% 
  summarise(Count = sum(Count)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = 'ID', values_from = 'Count', values_fill = list(Count = 0)) 

abundanceWide <- abundanceWide %>% 
  mutate(sum = rowSums(abundanceWide[, c(8:42)], na.rm = TRUE)) %>% 
  filter(sum > 0, !(Sample %in% c('C1A_A_38_ABG', 'C3SA_D_16_PBG', 'C3SA_C_38_PBG')))


print(permanova <- adonis2(formula = abundanceWide[,8:42]~Treatment, data=abundanceWide, permutations=999, method="bray"))
#F=2.3372   , df=1,30, p=0.016 * 

#betadisper
veg <- vegdist(abundanceWide[,8:42], method = "bray")
dispersion <- betadisper(veg, abundanceWide$Treatment)
permutest(dispersion, pairwise=TRUE, permutations=999) 

#F=2.0478   , df=1,30, p=0.159

BC_Data <- metaMDS(abundanceWide[,8:42])
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
ABG_VS_PBG_NMDS_2023 <- ggplot(BC_NMDS_Graph, aes(x = MDS1, y = MDS2, color = group, linetype = group, shape = group)) +
  geom_point(size = 6) + 
  geom_path(data = BC_Ellipses, aes(x = NMDS1, y = NMDS2), size = 3) +
  labs(color = "", linetype = "", shape = "") +
  scale_colour_manual(values = c("blue", "red"), name = "") +
  scale_linetype_manual(values = c("solid", "twodash"), name = "") +
  scale_shape_manual(values = c(19, 19)) +
  xlab("NMDS1") + 
  ylab("NMDS2") + 
  theme_classic() +
  theme(axis.text.x = element_text(size = NMDS_AXIS_TEXT_SIZE, color = "black"), 
        axis.text.y = element_text(size = NMDS_AXIS_TEXT_SIZE, color = "black"), 
        axis.title.x = element_text(size = NMDS_AXIS_TEXT_SIZE, color = 'black'),
        axis.title.y = element_text(size = NMDS_AXIS_TEXT_SIZE, color = 'black'),
        legend.text = element_text(size = 60),
        legend.key.size = unit(5, "cm"),  # Adjust the size of legend key (dots)
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        plot.title = element_text(size = 75)
  ) +
  ggtitle("2023 ABG vs PBG") + 
  annotate('text', x = min(BC_NMDS_Graph$MDS1) - 0.04, y = c(min(BC_NMDS_Graph$MDS2) + 0.02, min(BC_NMDS_Graph$MDS2) - 0.2), 
           label = c('Mean p = 0.016', 'Variance p = 0.159'), size = STATS_TEXT_SIZE, hjust = 'left')

#mean F=2.3372   , df=1,30, p=0.016 * 
#variance F=2.0478   , df=1,30, p=0.159

ABG_VS_PBG_NMDS_2023

#### Multipanel NMDS ####



NMMDS_BIG <- plot_grid(
  ABG_VS_PBG_NMDS_2021, Years_Since_Burned_NMDS_2021,
  ABG_VS_PBG_NMDS_2022, Years_Since_Burned_NMDS_2022,
  ABG_VS_PBG_NMDS_2023, Years_Since_Burned_NMDS_2023,
  ncol = 2
)


# NMMDS_BIG <- grid.arrange(
#   ABG_VS_PBG_NMDS_2021, Years_Since_Burned_NMDS_2021,
#   ABG_VS_PBG_NMDS_2022, Years_Since_Burned_NMDS_2022,
#   ABG_VS_PBG_NMDS_2023, Years_Since_Burned_NMDS_2023,
#   ncol = 2
# )


NMMDS_ABG_VS_PBG <- grid.arrange(
 ABG_VS_PBG_NMDS_2021, 
 ABG_VS_PBG_NMDS_2022, 
 ABG_VS_PBG_NMDS_2023,
  ncol = 1
)

ggsave("NMMDS_ABG_VS_PBG.png", NMMDS_ABG_VS_PBG, width = 40, height = 40)

ggsave("2021-2023 NMDS.png", NMMDS_BIG, width = 40, height = 40)

#### Overall PERMANOVA & NMDS ####

combined_data <- rbind(Abundance_Data2021, Abundance_Data2022, Abundance_Data2023)

combined_data_wide <- combined_data %>% 
  mutate(block = ifelse(grepl("S", Sample), "North", "South")) %>%
  select(Sample, block, WS, Trans, Plot, Treatment, TreatmentSB, Count, ID) %>%
  group_by(Sample, block, WS, Trans, Plot, Treatment, TreatmentSB, ID) %>% 
  summarise(Count = sum(Count)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = 'ID', values_from = 'Count', values_fill = list(Count = 0)) 

combined_data_wide <- combined_data_wide %>% 
  mutate(sum = rowSums(combined_data_wide[, c(8:124)], na.rm = TRUE)) %>% 
  filter(sum > 0, !(Sample %in% c('C1A_A_38_ABG', 'C3SA_D_16_PBG', 'C3SA_C_38_PBG')))

#Years since burned
print(permanova <- adonis2(
  combined_data_wide[,8:124]~TreatmentSB,
  data = combined_data_wide,
  method = "bray",
  permutations=999,
  strata = combined_data_wide$Sample))

#F = 2.3383, P = 0.127

#betadisper
veg <- vegdist(combined_data_wide[,8:124], method = "bray")
dispersion <- betadisper(veg, combined_data_wide$TreatmentSB)
permutest(dispersion, pairwise=TRUE, permutations=999) 
#F=0.4157 , df=3,191, p=0.739

BC_Data <- metaMDS(combined_data_wide[,8:124])
sites <- 1:nrow(combined_data_wide)
BC_Meta_Data <- combined_data_wide[,1:7]
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
  theme_bw() +
  theme(axis.text.x=element_text(size=24, color = "black"), 
        axis.text.y = element_text(size = 24, color = "black"), 
        axis.title.x = element_text(size = 24, color = 'black'),
        axis.title.y = element_text(size = 24, color = 'black'),
        legend.text = element_text(size = 24),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()
  ) +
  ggtitle("all years, year since burned")

#F=1.6053, df=3,52, p=0.016
#export at 1500x1000


#Add NMDS

#### Prep for Bootstrapping ####

#2021 Prep
Abundance2021 <- Abundance_ID_aboveground %>% filter(Date == 2021) %>% 
  filter(Plot %in% c("2", "4"))
commMetircs2021 <-  Abundance2021 %>% community_structure(abundance.var='Count', replicate.var='Sample')
commMetircs2021$Treatment <- ifelse(grepl("C1", commMetircs2021$Sample), "ABG", "PBG") 
commMetircs2021 <- commMetircs2021 %>% separate(Sample, into = c('WS', 'Trans', 'Plot'), sep = '_')

#making a 2021 total counts
total_counts_2021 <- aggregate(Count ~ Sample + Treatment, data = Abundance2021, FUN = sum) %>% 
  separate(Sample, into = c("WS", "Trans", "Plot"), sep = "_", extra = "drop") %>% 
  mutate(block = ifelse(grepl("S", WS), "North", "South"))

#2022 Prep
Abundance2022 <- Abundance_ID_aboveground %>% filter(Date == 2022)
commMetircs2022 <-  Abundance2022 %>% community_structure(abundance.var='Count', replicate.var='Sample')
commMetircs2022$Treatment <- ifelse(grepl("C1", commMetircs2022$Sample), "ABG", "PBG") 
commMetircs2022 <- commMetircs2022 %>% separate(Sample, into = c('WS', 'Trans', 'Plot'), sep = '_')

#making a 2022 total counts
total_counts_2022 <- aggregate(Count ~ Sample + Treatment, data = Abundance2022, FUN = sum) %>% 
  separate(Sample, into = c("WS", "Trans", "Plot"), sep = "_", extra = "drop") %>% 
  mutate(block = ifelse(grepl("S", WS), "North", "South"))

#2023 Prep
Abundance2023 <- Abundance_ID_aboveground %>% filter(Date == 2023)
commMetircs2023 <-  Abundance2023 %>% community_structure(abundance.var='Count', replicate.var='Sample')
commMetircs2023$Treatment <- ifelse(grepl("C1", commMetircs2023$Sample), "ABG", "PBG") 
commMetircs2023 <- commMetircs2023 %>% separate(Sample, into = c('WS', 'Trans', 'Plot'), sep = '_')

#making a 2023 total counts
total_counts_2023 <- aggregate(Count ~ Sample + Treatment, data = Abundance2023, FUN = sum) %>% 
  separate(Sample, into = c("WS", "Trans", "Plot"), sep = "_", extra = "drop") %>% 
  mutate(block = ifelse(grepl("S", WS), "North", "South"))


#Making Joined Dataframs for all years
Joined2021 <- full_join(BurnInfo2021, commMetircs2021) %>% unite("TreatmentSB",c("Treatment","SB"), sep="_") %>% 
  mutate(block = ifelse(grepl("S", WS), "North", "South")) %>% 
  mutate(Treatment = ifelse(grepl("C1", WS), "ABG", "PBG"))

Joined2022 <- full_join(BurnInfo2022, commMetircs2022) %>% unite("TreatmentSB",c("Treatment","SB"), sep="_") %>% 
  mutate(block = ifelse(grepl("S", WS), "North", "South")) %>% 
  mutate(Treatment = ifelse(grepl("C1", WS), "ABG", "PBG"))

Joined2023 <- full_join(BurnInfo2023, commMetircs2023) %>% unite("TreatmentSB",c("Treatment","SB"), sep="_") %>% 
  mutate(block = ifelse(grepl("S", WS), "North", "South")) %>% 
  mutate(Treatment = ifelse(grepl("C1", WS), "ABG", "PBG"))

#### 2021 Bootstrapping! ####
set.seed(123)

Weight_2021 <- Weight_2_combined %>%
  separate(sample, into = c("Date", "WS", "Trans", "Plot"), sep = "_") %>%
  group_by(block) %>% 
  filter(Date == "2021") %>% 
  filter(Treatment == "PBG") %>% 
  filter(Plot %in% c("2", "4")) %>% 
  mutate(plot_index=1:length(block)) %>% 
  select(-Date)


Joined_New_2021 <- Joined2021 %>%
  filter(Treatment == "PBG") %>%
  group_by(block) %>%
  mutate(plot_index=1:length(block)) %>% 
  filter(Plot %in% c("2", "4"))

total_counts_2021 <- total_counts_2021 %>% filter(Treatment == "PBG") %>% group_by(block) %>% 
  mutate(plot_index=1:length(block)) %>%  filter() %>%  
  filter(Plot %in% c("2", "4"))

Combined_Data_2021 <- left_join(Joined_New_2021, total_counts_2021, by = c("WS", "Trans", "Plot", "block", "plot_index")) %>% 
  select(-TreatmentSB) %>%
  rename(Treatment = Treatment.x) %>%
  select(-Treatment.y)

Combined_Data_2021 <- left_join(Weight_2021, Combined_Data_2021, by = c("WS", "Trans", "Plot", "block", "plot_index"))


#Drop big outliers 
Combined_Data_2021 <- Combined_Data_2021 %>%
  filter(combined_weight <= 190)


num_bootstrap <- 1000
bootstrap_vector <- 1:num_bootstrap
PBG_plot_master_2021 <- data.frame()  # Initialize an empty dataframe

for (BOOT in bootstrap_vector) {
  Joined_New_Key_2021 <- Combined_Data_2021 %>%
    dplyr::select(1:10) %>%
    unique() %>%
    group_by(block) %>%
    sample_n(16, replace = TRUE) %>%
    dplyr::select(plot_index) %>%
    ungroup()
  
  # Join the sampled rows back to the original dataframe
  PBG_plot_ready_2021 <- Combined_Data_2021 %>%
    right_join(Joined_New_Key_2021, by = c("block", "plot_index")) %>%
    mutate(iteration = BOOT)
  
  # Append the results to the master dataframe
  PBG_plot_master_2021 <- rbind(PBG_plot_master_2021, PBG_plot_ready_2021)
}

#### 2022 Bootstrapping! ####
set.seed(123)


Weight_2022 <- Weight_2_combined %>%
  separate(sample, into = c("Date", "WS", "Trans", "Plot"), sep = "_") %>%
  group_by(block) %>% 
  filter(Date == "2022") %>% 
  filter(Treatment == "PBG") %>% 
  select(-Date)

Joined_New_2022 <- Joined2022 %>%
  filter(Treatment == "PBG") %>%
  group_by(block) %>%
  mutate(plot_index=1:length(block))


total_counts_2022 <- total_counts_2022 %>% filter(Treatment == "PBG") %>% 
  group_by(block) %>% 
  mutate(plot_index=1:length(block)) %>%  filter()



Combined_Data_2022 <- left_join(Joined_New_2022, total_counts_2022, by = c("WS", "Trans", "Plot", "block", "plot_index")) %>% 
  select(-TreatmentSB) %>%
  rename(Treatment = Treatment.x) %>%
  select(-Treatment.y)




Combined_Data_2022 <- left_join(Weight_2022, Combined_Data_2022, by = c("WS", "Trans", "Treatment", "Plot", "block"))

Combined_Data_2022 <- na.omit(Combined_Data_2022)

num_bootstrap <- 1000
bootstrap_vector <- 1:num_bootstrap
PBG_plot_master_2022<- data.frame()  # Initialize an empty dataframe

for (BOOT in bootstrap_vector) {
  Joined_New_Key_2022 <- Combined_Data_2022 %>%
    dplyr::select(1:10) %>%
    unique() %>%
    group_by(block) %>%
    sample_n(16, replace = TRUE) %>%
    dplyr::select(plot_index) %>%
    ungroup()
  
  # Join the sampled rows back to the original dataframe
  PBG_plot_ready_2022 <- Combined_Data_2022 %>%
    right_join(Joined_New_Key_2022, by = c("block", "plot_index")) %>%
    mutate(iteration = BOOT)
  
  # Append the results to the master dataframe
  PBG_plot_master_2022 <- rbind(PBG_plot_master_2022, PBG_plot_ready_2022)
}

#### 2023 Bootstrapping! ####
set.seed(123)

Weight_2023 <- Weight_2_combined %>%
  separate(sample, into = c("Date", "WS", "Trans", "Plot"), sep = "_") %>%
  group_by(block) %>% 
  filter(Date == "2023") %>% 
  filter(Treatment == "PBG") %>% 
  filter(Plot %in% c("1", "3")) %>% 
  mutate(plot_index=1:length(block)) %>% 
  select(-Date)


Joined_New_2023 <- Joined2023 %>%
  filter(Treatment == "PBG") %>%
  group_by(block) %>%
  mutate(plot_index=1:length(block)) %>% 
  filter(Plot %in% c("1", "3"))

total_counts_2023 <- total_counts_2023 %>% filter(Treatment == "PBG") %>% group_by(block) %>% 
  mutate(plot_index=1:length(block)) %>%  filter() %>%  
  filter(Plot %in% c("1", "3"))

Combined_Data_2023 <- left_join(Joined_New_2023, total_counts_2023, by = c("WS", "Trans", "Plot", "block", "plot_index")) %>% 
  select(-TreatmentSB) %>%
  rename(Treatment = Treatment.x) %>%
  select(-Treatment.y)

Combined_Data_2023 <- left_join(Weight_2023, Combined_Data_2023, by = c("WS", "Trans", "Plot", "block", "plot_index"))



num_bootstrap <- 1000
bootstrap_vector <- 1:num_bootstrap
PBG_plot_master_2023 <- data.frame()  # Initialize an empty dataframe

for (BOOT in bootstrap_vector) {
  Joined_New_Key_2023 <- Combined_Data_2023 %>%
    dplyr::select(1:10) %>%
    unique() %>%
    group_by(block) %>%
    sample_n(16, replace = TRUE) %>%
    dplyr::select(plot_index) %>%
    ungroup()
  
  # Join the sampled rows back to the original dataframe
  PBG_plot_ready_2023 <- Combined_Data_2023 %>%
    right_join(Joined_New_Key_2023, by = c("block", "plot_index")) %>%
    mutate(iteration = BOOT)
  
  # Append the results to the master dataframe
  PBG_plot_master_2023 <- rbind(PBG_plot_master_2023, PBG_plot_ready_2023)
}


#### 2021 Z-Score  Calculations ####


#Getting average richness per iteration for bootstrapped dataframe
average_richness_2021 <- PBG_plot_master_2021 %>%
  group_by(iteration) %>%
  summarize(mean_richness = mean(richness))

#Getting evenness richness per iteration for bootstrapped dataframe
average_evenness_2021 <- PBG_plot_master_2021 %>%
  group_by(iteration) %>%
  summarize(mean_evenness = mean(Evar, na.rm = TRUE))

#Getting average richness per iteration for bootstrapped dataframe
average_total_count_2021 <- PBG_plot_master_2021 %>%
  group_by(iteration) %>%
  summarize(mean_count = mean(Count, na.rm = TRUE))

#Getting average weight per iteration for bootstrapped dataframe
average_total_weight_2021 <- PBG_plot_master_2021 %>%
  group_by(iteration) %>%
  summarize(mean_weight = mean(combined_weight, na.rm = TRUE))

#Take the mean of the mean for PBG richness
PBG_mean_mean_richness_2021 <- mean(average_richness_2021$mean_richness, na.rm = TRUE)

#Take the mean of the mean for PBG evenness
PBG_mean_mean_evenness_2021 <- mean(average_evenness_2021$mean_evenness, na.rm = TRUE)

#Take the mean of the mean for PBG total count
PBG_mean_mean_total_count_2021 <- mean(average_total_count_2021$mean_count, na.rm = TRUE)

#Take the mean of the mean for PBG weight
PBG_mean_mean_total_weight_2021 <- mean(average_total_weight_2021$mean_weight, na.rm = TRUE)

#Take the mean of the mean for PBG total count

#Getting ABG ready

### 2021 ABG Count ###

Z_total_counts_2021 <- Abundance_ID_aboveground %>%
  filter(Date == 2021) %>%
  group_by(WS, Trans, Plot, Treatment) %>%
  summarise(Count = sum(Count, na.rm = TRUE)) %>% 
  filter(Plot %in% c("2", "4"))

total_counts_ABG_2021 <- Z_total_counts_2021 %>% filter(Treatment == "ABG")

#Take the mean of the mean for ABG total count
ABG_mean_total_count_2021 <- mean(total_counts_ABG_2021$Count, na.rm = TRUE)

### 2021 ABG Richness ###

Z_total_richness_2021 <- commMetircs2021 %>%
  group_by(WS, Trans, Plot, Treatment) %>%
  summarise(richness = sum(richness, na.rm = TRUE))

Total_richness_ABG_2021 <- Z_total_richness_2021 %>% filter(Treatment == "ABG")

#Take the mean of the mean for ABG total count
mean_richness_ABG_2021 <- mean(Total_richness_ABG_2021$richness, na.rm = TRUE)

### 2021 ABG Evenness ###

Z_total_evenness_2021 <- commMetircs2021 %>%
  group_by(WS, Trans, Plot, Treatment) %>%
  summarise(Evar = sum(Evar, na.rm = TRUE))

Total_evenness_ABG_2021 <- Z_total_evenness_2021 %>% filter(Treatment == "ABG")

#Take the mean of the mean for ABG total count
mean_evenness_ABG_2021 <- mean(Total_evenness_ABG_2021$Evar, na.rm = TRUE)

### 2021 ABG Weight
Z_total_weight_2021 <- Weight_2_combined %>%  
  separate(sample, into = c("Date", "WS", "Trans", "Plot"), sep = "_") %>% 
  group_by(WS, Trans, Plot, Treatment) %>%
  filter(Date == "2021") %>% 
  summarise(combined_weight = sum(combined_weight, na.rm = TRUE))

Weight_ABG_2021 <- Z_total_weight_2021 %>% filter(Treatment == "ABG")

Weight_ABG_2021 <- Weight_ABG_2021 %>%
  filter(combined_weight <= 190)

#Take the mean of the mean for ABG weight
mean_weight_ABG_2021 <- mean(Weight_ABG_2021$combined_weight, na.rm = TRUE)

### Z-Score for richness 

Z_R_2021 <- ((mean_richness_ABG_2021) - (PBG_mean_mean_richness_2021))/(sd(average_richness_2021$mean_richness))
Z_R_2021

p_value_R_2021 <- 2*pnorm(-abs(Z_R_2021))

print(p_value_R_2021)

#1.235378
#P = 0.2166899

# Z-Score for evenness 
Z_E_2021 <- ((mean_evenness_ABG_2021) - (PBG_mean_mean_evenness_2021))/(sd(average_evenness_2021$mean_evenness))
Z_E_2021

p_value_E_2021 <- 2*pnorm(-abs(Z_E_2021))

print(p_value_E_2021)
#--0.1450691
#P = 0.8846563

# Z-Score for total count 
Z_C_2021 <- ((ABG_mean_total_count_2021) - (PBG_mean_mean_total_count_2021))/(sd(average_total_count_2021$mean_count))
Z_C_2021

p_value_C_2021 <- 2*pnorm(-abs(Z_C_2021))

print(p_value_C_2021)
#Z-Score: 1.145737
#P = 0.2519041

# Z-Score for Weight

Z_W_2021 <- ((mean_weight_ABG_2021) - (PBG_mean_mean_total_weight_2021))/(sd(average_total_weight_2021$mean_weight))
Z_W_2021

p_value_W_2021 <- 2*pnorm(-abs(Z_W_2021))

print(p_value_W_2021)

#Z-score = -0.8275585
# P = 0.4079206

#### 2022 Z-Score Calculations ####

#Getting average richness per iteration for bootstrapped dataframe
average_richness_2022 <- PBG_plot_master_2022 %>%
  group_by(iteration) %>%
  summarize(mean_richness = mean(richness))

#Getting evenness richness per iteration for bootstrapped dataframe
average_evenness_2022 <- PBG_plot_master_2022 %>%
  group_by(iteration) %>%
  summarize(mean_evenness = mean(Evar, na.rm = TRUE))

#Getting average richness per iteration for bootstrapped dataframe
average_total_count_2022 <- PBG_plot_master_2022 %>%
  group_by(iteration) %>%
  summarize(mean_count = mean(Count, na.rm = TRUE))

#Getting average weight per iteration for bootstrapped dataframe
average_total_weight_2022<- PBG_plot_master_2022 %>%
  group_by(iteration) %>%
  summarize(mean_weight = mean(combined_weight, na.rm = TRUE))

#Take the mean of the mean for PBG richness
PBG_mean_mean_richness_2022 <- mean(average_richness_2022$mean_richness, na.rm = TRUE)

#Take the mean of the mean for PBG evenness
PBG_mean_mean_evenness_2022 <- mean(average_evenness_2022$mean_evenness, na.rm = TRUE)

#Take the mean of the mean for PBG evenness
PBG_mean_mean_total_count_2022 <- mean(average_total_count_2022$mean_count, na.rm = TRUE)

#Take the mean of the mean for PBG weight
PBG_mean_mean_total_weight_2022 <- mean(average_total_weight_2022$mean_weight, na.rm = TRUE)


#Getting ABG ready

### 2022 ABG Count ###

Z_total_counts_2022 <- Abundance_ID_aboveground %>%
  filter(Date == 2022) %>%
  group_by(WS, Trans, Plot, Treatment) %>%
  summarise(Count = sum(Count, na.rm = TRUE))

total_counts_ABG_2022 <- Z_total_counts_2022 %>% filter(Treatment == "ABG")

#Take the mean of the mean for ABG total count
ABG_mean_total_count_2022 <- mean(total_counts_ABG_2022$Count, na.rm = TRUE)

### 2022 ABG Richness ###

Z_total_richness_2022 <- commMetircs2022 %>%
  group_by(WS, Trans, Plot, Treatment) %>%
  summarise(richness = sum(richness, na.rm = TRUE))

Total_richness_ABG_2022 <- Z_total_richness_2022 %>% filter(Treatment == "ABG")

#Take the mean of the mean for ABG total count
mean_richness_ABG_2022 <- mean(Total_richness_ABG_2022$richness, na.rm = TRUE)

### 2022 ABG Evenness ###

Z_total_evenness_2022 <- commMetircs2022 %>%
  group_by(WS, Trans, Plot, Treatment) %>%
  summarise(Evar = sum(Evar, na.rm = TRUE))

Total_evenness_ABG_2022 <- Z_total_evenness_2022 %>% filter(Treatment == "ABG")

#Take the mean of the mean for ABG total count
mean_evenness_ABG_2022 <- mean(Total_evenness_ABG_2022$Evar, na.rm = TRUE)

### 2021 ABG Weight
Z_total_weight_2022 <- Weight_2_combined %>%  
  separate(sample, into = c("Date", "WS", "Trans", "Plot"), sep = "_") %>% 
  group_by(WS, Trans, Plot, Treatment) %>%
  filter(Date == "2022") %>% 
  summarise(combined_weight = sum(combined_weight, na.rm = TRUE))

Weight_ABG_2022 <- Z_total_weight_2022 %>% filter(Treatment == "ABG")

#Take the mean of the mean for ABG weight
mean_weight_ABG_2022 <- mean(Weight_ABG_2022$combined_weight, na.rm = TRUE)

# Z-Score for richness 
Z_R_2022 <- ((mean_richness_ABG_2022) - (PBG_mean_mean_richness_2022))/(sd(average_richness_2022$mean_richness))
Z_R_2022
p_value_R_2022 <- 2*pnorm(-abs(Z_R_2022))
print(p_value_R_2022)
#2.028526
#P = 0.04250654

# Z-Score for evenness 
Z_E_2022 <- ((mean_evenness_ABG_2022) - (PBG_mean_mean_evenness_2022))/(sd(average_evenness_2022$mean_evenness))
Z_E_2022
p_value_E_2022 <- 2*pnorm(-abs(Z_E_2022))
print(p_value_E_2022)
#-0.9478015
#P = 0.3432305

# Z-Score for total count 
Z_C_2022 <- ((ABG_mean_total_count_2022) - (PBG_mean_mean_total_count_2022))/(sd(average_total_count_2022$mean_count))
Z_C_2022
p_value_C_2022 <- 2*pnorm(-abs(Z_C_2022))
print(p_value_C_2022)
#Z-Score: 1.514056
#P = 0.1300116

#Z-score for Weight
Z_W_2022 <- ((mean_weight_ABG_2022) - (PBG_mean_mean_total_weight_2022))/(sd(average_total_weight_2022$mean_weight))
Z_W_2022

p_value_W_2022 <- 2*pnorm(-abs(Z_W_2022))

print(p_value_W_2022)

#Z-Score: 1.69659
#P-Value: 0.08977415

#### 2023 Z-Score Calculations ####
#Getting average richness per iteration for bootstrapped dataframe
average_richness_2023 <- PBG_plot_master_2023 %>%
  group_by(iteration) %>%
  summarize(mean_richness = mean(richness))

#Getting evenness richness per iteration for bootstrapped dataframe
average_evenness_2023 <- PBG_plot_master_2023 %>%
  group_by(iteration) %>%
  summarize(mean_evenness = mean(Evar, na.rm = TRUE))

#Getting average richness per iteration for bootstrapped dataframe
average_total_count_2023 <- PBG_plot_master_2023 %>%
  group_by(iteration) %>%
  summarize(mean_count = mean(Count, na.rm = TRUE))

#Getting average weight per iteration for bootstrapped dataframe
average_total_weight_2023 <- PBG_plot_master_2023 %>%
  group_by(iteration) %>%
  summarize(mean_weight = mean(combined_weight, na.rm = TRUE))

#Take the mean of the mean for PBG richness
PBG_mean_mean_richness_2023 <- mean(average_richness_2023$mean_richness, na.rm = TRUE)

#Take the mean of the mean for PBG evenness
PBG_mean_mean_evenness_2023 <- mean(average_evenness_2023$mean_evenness, na.rm = TRUE)

#Take the mean of the mean for PBG evenness
PBG_mean_mean_total_count_2023 <- mean(average_total_count_2023$mean_count, na.rm = TRUE)

#Take the mean of the mean for PBG weight
PBG_mean_mean_total_weight_2023 <- mean(average_total_weight_2023$mean_weight, na.rm = TRUE)

#Getting ABG ready

### 2023 ABG Count ###

Z_total_counts_2023 <- Abundance_ID_aboveground %>%
  filter(Date == 2023) %>%
  group_by(WS, Trans, Plot, Treatment) %>%
  summarise(Count = sum(Count, na.rm = TRUE))

total_counts_ABG_2023 <- Z_total_counts_2023 %>% filter(Treatment == "ABG")

#Take the mean of the mean for ABG total count
ABG_mean_total_count_2023 <- mean(total_counts_ABG_2023$Count, na.rm = TRUE)

### 2023 ABG Richness ###

Z_total_richness_2023 <- commMetircs2023 %>%
  group_by(WS, Trans, Plot, Treatment) %>%
  summarise(richness = sum(richness, na.rm = TRUE))

Total_richness_ABG_2023 <- Z_total_richness_2023 %>% filter(Treatment == "ABG")

#Take the mean of the mean for ABG total count
mean_richness_ABG_2023 <- mean(Total_richness_ABG_2023$richness, na.rm = TRUE)

### 2023 ABG Evenness ###

Z_total_evenness_2023 <- commMetircs2023 %>%
  group_by(WS, Trans, Plot, Treatment) %>%
  summarise(Evar = sum(Evar, na.rm = TRUE))

Total_evenness_ABG_2023 <- Z_total_evenness_2023 %>% filter(Treatment == "ABG")

### 2021 ABG Weight
Z_total_weight_2023 <- Weight_2_combined %>%  
  separate(sample, into = c("Date", "WS", "Trans", "Plot"), sep = "_") %>% 
  group_by(WS, Trans, Plot, Treatment) %>%
  filter(Date == "2023") %>% 
  summarise(combined_weight = sum(combined_weight, na.rm = TRUE))

Weight_ABG_2023 <- Z_total_weight_2023 %>% filter(Treatment == "ABG")

#Take the mean of the mean for ABG weight
mean_weight_ABG_2023 <- mean(Weight_ABG_2023$combined_weight, na.rm = TRUE)

#Take the mean of the mean for ABG total count
mean_evenness_ABG_2023 <- mean(Total_evenness_ABG_2023$Evar, na.rm = TRUE)

# Z-Score for richness 
Z_R_2023 <- ((mean_richness_ABG_2023) - (PBG_mean_mean_richness_2023))/(sd(average_richness_2023$mean_richness))
Z_R_2023
p_value_R_2023 <- 2*pnorm(-abs(Z_R_2023))
print(p_value_R_2023)
#-0.4446148
#P =  0.6565981

# Z-Score for evenness 
Z_E_2023 <- ((mean_evenness_ABG_2023) - (PBG_mean_mean_evenness_2023))/(sd(average_evenness_2023$mean_evenness))
Z_E_2023
p_value_E_2023 <- 2*pnorm(-abs(Z_E_2023))
print(p_value_E_2023)
#-3.400096
#P =0.000673623

# Z-Score for total count 
Z_C_2023 <- ((ABG_mean_total_count_2023) - (PBG_mean_mean_total_count_2023))/(sd(average_total_count_2023$mean_count))
Z_C_2023
p_value_C_2023 <- 2*pnorm(-abs(Z_C_2023))
print(p_value_C_2023)
#Z-Score: 0.1427381
#P = 0.8864971

#Z-score for Weight
Z_W_2023<- ((mean_weight_ABG_2023) - (PBG_mean_mean_total_weight_2023))/(sd(average_total_weight_2023$mean_weight))
Z_W_2023

p_value_W_2023 <- 2*pnorm(-abs(Z_W_2023))

print(p_value_W_2023)

#Z-Score = 1.106847
#P-value = 0.2683601


#### Graphs Bootstrapped Means 2021 ####

#Richness means graph
richness_2021 <- ggplot(average_richness_2021, aes(x = mean_richness, y = ..scaled.., color = "PBG")) +
  geom_density(alpha = 0.5) +
  geom_vline(aes(xintercept = mean_richness_ABG_2021, color = "ABG"), linetype = "dashed", size = 1) +
  labs(title = "2021: Density Plot of Mean Richness",
       x = "Mean Richness",
       y = "Density") +
  scale_color_manual(values = c("blue", "red"), name = "Legend") +
  scale_x_continuous(limits = c(5, 20)) +
  scale_y_continuous(limits = c(0, 1)) +
  theme(legend.position = c(0.15, 0.7),
        axis.text.x = element_text(size = 22),  # Adjust the size as needed
        axis.text.y = element_text(size = 22))  # Adjust the size as needed

richness_2021



#Evenness means graph
evenness_2021 <- ggplot(average_evenness_2021, aes(x = mean_evenness, y = ..scaled.., color = "PBG")) +
  geom_density(alpha = 0.5) +
  geom_vline(aes(xintercept = mean_evenness_ABG_2021, color = "ABG Mean Evenness"), linetype = "dashed", size = 1) +
  labs(title = "2021: Density Plot of Mean Evenness",
       x = "Mean Evenness",
       y = "Density") +
  scale_color_manual(values = c("blue", "red"), name = "Legend") +
#  annotate("text", x = 0.8, y = 1, label = "p-value: 0.885", color = "blue", size = 4, hjust = 1, vjust = 1) +
#  annotate("text", x = 0.8, y = 0.8, label = "Z-score: 0.145", color = "red", size = 4, hjust = 1, vjust = 1) +
  theme(legend.position = "none",  plot.title = element_text(size = 15), 
        axis.text.x = element_text(size = 22),  # Adjust the size as needed
        axis.text.y = element_text(size = 22)) +
  scale_x_continuous(limits = c(.5, .8)) +
  scale_y_continuous(limits = c(0, 1)) 

evenness_2021

#Total_count means graphs
count_2021 <- ggplot(average_total_count_2021, aes(x = mean_count, y = ..scaled.., color = "PBG")) +
  geom_density(alpha = 0.5) +
  geom_vline(aes(xintercept = ABG_mean_total_count_2021, color = "ABG Total Abundance"), linetype = "dashed", size = 1) +
  labs(title = "2021: Density Plot of Total Abundance",
       x = "Mean Total Abundance",
       y = "Density") +
  scale_color_manual(values = c("blue", "red"), name = "Legend") +
#  annotate("text", x = 55, y = 1, label = "p-value: 0.252", color = "blue", size = 4, hjust = 1, vjust = 1) +
#  annotate("text", x = 55, y = 0.80, label = "Z-score: 1.15", color = "red", size = 4, hjust = 1, vjust = 1) +
  theme(legend.position = "none",  plot.title = element_text(size = 15),  axis.text.x = element_text(size = 22),  # Adjust the size as needed
        axis.text.y = element_text(size = 22)) +
  scale_x_continuous(limits = c(15, 55)) +
  scale_y_continuous(limits = c(0, 1)) 

count_2021

#Total_weight means graphs
weight_2021 <- ggplot(average_total_weight_2021, aes(x = mean_weight, y = ..scaled.., color = "PBG")) +
  geom_density(alpha = 0.5) +
  geom_vline(aes(xintercept = mean_weight_ABG_2021, color = "ABG Total Count"), linetype = "dashed", size = 1) +
  labs(title = "2021: Density Plot of Weight",
       x = "Mean Weight",
       y = "Density") +
  scale_color_manual(values = c("blue", "red"), name = "Legend") +
#  annotate("text", x = 130, y = 1, label = "p-value: 0.408", color = "blue", size = 4, hjust = 1, vjust = 1) +
#  annotate("text", x = 130, y = 0.8, label = "Z-score: 0.828", color = "red", size = 4, hjust = 1, vjust = 1) +
  theme(legend.position = "none",  plot.title = element_text(size = 15),  axis.text.x = element_text(size = 22),  # Adjust the size as needed
        axis.text.y = element_text(size = 22)) +
  scale_x_continuous(limits = c(25, 130)) +
  scale_y_continuous(limits = c(0, 1)) 

weight_2021

grid.arrange(richness_2021, evenness_2021, count_2021, weight_2021)

#### Graphs Bootstrapped Means 2022 ####
#Richness means graph
richness_2022 <- ggplot(average_richness_2022, aes(x = mean_richness, y = ..scaled..,  color = "PBG")) +
  geom_density(alpha = 0.5) +
  geom_vline(aes(xintercept = mean_richness_ABG_2022 , color = "ABG Mean Richness"), linetype = "dashed", size = 1) +
  labs(title = "2022: Density Plot of Mean Richness",
       x = "Mean Richness",
       y = "Density") +
  scale_color_manual(values = c("blue", "red"), name = "Legend")  +
#  annotate("text", x = 20, y = 1, label = "p-value: 0.0425", color = "blue", size = 4, hjust = 1, vjust = 1) +
#  annotate("text", x = 20, y = .80, label = "Z-score: 2.03", color = "red", size = 4, hjust = 1, vjust = 1) +
  theme(legend.position = "none",  plot.title = element_text(size = 15), 
        axis.text.x = element_text(size = 22),  # Adjust the size as needed
        axis.text.y = element_text(size = 22)) +
  scale_x_continuous(limits = c(5, 20)) +
  scale_y_continuous(limits = c(0, 1)) 
richness_2022

#Evenness means graph
evenness_2022 <- ggplot(average_evenness_2022, aes(x = mean_evenness, y = ..scaled..,  color = "PBG")) +
  geom_density(alpha = 0.5) +
  geom_vline(aes(xintercept = mean_evenness_ABG_2022, color = "ABG Mean Evenness"), linetype = "dashed", size = 1) +
  labs(title = "2022: Density Plot of Mean Evenness",
       x = "Mean Evenness",
       y = "Density") +
  scale_color_manual(values = c("blue", "red"), name = "Legend") +
#  annotate("text", x = 0.8, y = 1, label = "p-value: 0.343", color = "blue", size = 4, hjust = 1, vjust = 1) +
#  annotate("text", x = 0.8, y = 0.8, label = "Z-score: 0.948", color = "red", size = 4, hjust = 1, vjust = 1) +
  theme(legend.position = "none",  plot.title = element_text(size = 15),axis.text.x = element_text(size = 22),  # Adjust the size as needed
        axis.text.y = element_text(size = 22) ) +
  scale_x_continuous(limits = c(.5, .8)) +
  scale_y_continuous(limits = c(0, 1)) 
evenness_2022

#Total_count means graphs
count_2022 <- ggplot(average_total_count_2022, aes(x = mean_count, y = ..scaled..,  color = "PBG")) +
  geom_density(alpha = 0.5) +
  geom_vline(aes(xintercept = ABG_mean_total_count_2022, color = "ABG Total Abundance"), linetype = "dashed", size = 1) +
  labs(title = "2022: Density Plot of Total Abundance",
       x = "Mean Total Abundance",
       y = "Density") +
  scale_color_manual(values = c("blue", "red"), name = "Legend") +
#  annotate("text", x = 55, y = 1, label = "p-value: 0.130", color = "blue", size = 4, hjust = 1, vjust = 1) +
#  annotate("text", x = 55, y = .80, label = "Z-score: 1.51", color = "red", size = 4, hjust = 1, vjust = 1) +
  theme(legend.position = "none",  plot.title = element_text(size = 15), axis.text.x = element_text(size = 22),  # Adjust the size as needed
        axis.text.y = element_text(size = 22) ) +
  scale_x_continuous(limits = c(15, 55)) +
  scale_y_continuous(limits = c(0, 1)) 
count_2022

#Total_weight means graphs
weight_2022 <- ggplot(average_total_weight_2022, aes(x = mean_weight, y = ..scaled..,  color = "PBG")) +
  geom_density(alpha = 0.5) +
  geom_vline(aes(xintercept = mean_weight_ABG_2022, color = "ABG Total Count"), linetype = "dashed", size = 1) +
  labs(title = "2022: Density Plot of Weight",
       x = "Mean Weight",
       y = "Density") +
  scale_color_manual(values = c("blue", "red"), name = "Legend") +
#  annotate("text", x = 130, y = 1, label = "p-value: 0.0897", color = "blue", size = 4, hjust = 1, vjust = 1) +
#  annotate("text", x = 130, y = .80, label = "Z-score: 1.697", color = "red", size = 4, hjust = 1, vjust = 1) +
  theme(legend.position = "none",  plot.title = element_text(size = 15),  axis.text.x = element_text(size = 22),  # Adjust the size as needed
        axis.text.y = element_text(size = 22)) +
  scale_x_continuous(limits = c(25, 130)) +
  scale_y_continuous(limits = c(0, 1)) 

weight_2022

grid.arrange(richness_2022, evenness_2022, count_2022)

#### Graphs Bootstrapped Means 2023 ####

#Richness means graph
richness_2023 <- ggplot(average_richness_2023, aes(x = mean_richness, y = ..scaled..,  color = "PBG")) +
  geom_density(alpha = 0.5) +
  geom_vline(aes(xintercept = mean_richness_ABG_2023 , color = "ABG Mean Richness"), linetype = "dashed", size = 1) +
  labs(title = "2023: Density Plot of Mean Richness",
       x = "Mean Richness",
       y = "Density") +
  scale_color_manual(values = c("blue", "red"), name = "Legend") +
#  annotate("text", x = 20, y = 1, label = "p-value: 0.657", color = "blue", size = 4, hjust = 1, vjust = 1) +
#  annotate("text", x = 20, y = 0.8, label = "Z-score: 0.445", color = "red", size = 4, hjust = 1, vjust = 1) +
  theme(legend.position = "none",  plot.title = element_text(size = 15),  axis.text.x = element_text(size = 22),  # Adjust the size as needed
        axis.text.y = element_text(size = 22)) +
  scale_x_continuous(limits = c(5, 20)) +
  scale_y_continuous(limits = c(0, 1)) 

#Evenness means graph
evenness_2023 <- ggplot(average_evenness_2023, aes(x = mean_evenness, y = ..scaled.., color = "PBG")) +
  geom_density(alpha = 0.5) +
  geom_vline(aes(xintercept = mean_evenness_ABG_2023, color = "ABG Mean Evenness"), linetype = "dashed", size = 1) +
  labs(title = "2023: Density Plot of Mean Evenness",
       x = "Mean Evenness",
       y = "Density") +
  scale_color_manual(values = c("blue", "red"), name = "Legend") +
#  annotate("text", x = 0.8, y = 1, label = "p-value: <0.001", color = "blue", size = 4, hjust = 1, vjust = 1) +
#  annotate("text", x = 0.8, y = 0.8, label = "Z-score:  3.400", color = "red", size = 4, hjust = 1, vjust = 1) +
  theme(legend.position = "none",  plot.title = element_text(size = 15),  axis.text.x = element_text(size = 22),  # Adjust the size as needed
        axis.text.y = element_text(size = 22)) +
  scale_x_continuous(limits = c(.5, .8)) +
  scale_y_continuous(limits = c(0, 1)) 



#Total_count means graphs
count_2023 <- ggplot(average_total_count_2023, aes(x = mean_count, y = ..scaled..,  color = "PBG")) +
  geom_density(alpha = 0.5) +
  geom_vline(aes(xintercept = ABG_mean_total_count_2023, color = "ABG Total Abundance"), linetype = "dashed", size = 1) +
  labs(title = "2023: Density Plot of Total Abundance",
       x = "Mean Total Abundance",
       y = "Density") +
  scale_color_manual(values = c("blue", "red"), name = "Legend") +
#  annotate("text", x = 55, y = 1, label = "p-value: 0.886", color = "blue", size = 4, hjust = 1, vjust = 1) +
#  annotate("text", x = 55, y = 0.8, label = "Z-score: 0.143", color = "red", size = 4, hjust = 1, vjust = 1) +
  theme(legend.position = "none",  plot.title = element_text(size = 15),  axis.text.x = element_text(size = 22),  # Adjust the size as needed
        axis.text.y = element_text(size = 22)) +
  scale_x_continuous(limits = c(15, 55)) +
  scale_y_continuous(limits = c(0, 1)) 

weight_2023 <- ggplot(average_total_weight_2023, aes(x = mean_weight, y = ..scaled..,  color = "PBG")) +
  geom_density(alpha = 0.5) +
  geom_vline(aes(xintercept = mean_weight_ABG_2023, color = "ABG Total Count"), linetype = "dashed", size = 1) +
  labs(title = "2023: Density Plot of Weight",
       x = "Mean Weight",
       y = "Density") +
  scale_color_manual(values = c("blue", "red")) +
#  annotate("text", x = 130, y = 1, label = "p-value: 0.268", color = "blue", size = 4, hjust = 1, vjust = 1) +
#  annotate("text", x = 130, y = 0.8, label = "Z-score: 1.11", color = "red", size = 4, hjust = 1, vjust = 1) +
  theme(legend.position = "none",  plot.title = element_text(size = 15),  axis.text.x = element_text(size = 22),  # Adjust the size as needed
        axis.text.y = element_text(size = 22)) +
  scale_x_continuous(limits = c(25, 130)) +
  scale_y_continuous(limits = c(0, 1)) 

weight_2023



grid.arrange(richness_2023, evenness_2023, count_2023)

#### Sample Graph For Legend ####
legend_2021 <- ggplot(average_richness_2021, aes(x = mean_richness, color = "PBG")) +
  geom_density(alpha = 0.5) +
  geom_vline(aes(xintercept = mean_richness_ABG_2021 , color = "ABG"), linetype = "dashed", size = 1) +
  labs(title = "2021: Density Plot of Mean Richness",
       x = "Mean Richness",
       y = "Density") +
  scale_color_manual(values = c("blue", "red"), name = "Legend") +
  annotate("text", x = 13, y = 0.5, label = "p-value: 0.202", color = "blue", size = 3) +
  annotate("text", x = 13, y = 0.44, label = "Z-score: 1.276", color = "red", size = 3) +
  theme_bw()

#### Big Graph of Means####

grid.arrange(
  richness_2021, evenness_2021, count_2021, weight_2021,
  richness_2022, evenness_2022, count_2022, weight_2022,
  richness_2023, evenness_2023, count_2023, weight_2023,
  ncol = 4,
  top = "Comparison of Richness, Evenness, and Count across Years"
)


# Arrange each plot individually
plot1 <- arrangeGrob(richness_2021, evenness_2021, count_2021, weight_2021, ncol = 1)
plot2 <- arrangeGrob(richness_2022, evenness_2022, count_2022, weight_2022, ncol = 1)
plot3 <- arrangeGrob(richness_2023, evenness_2023, count_2023, weight_2023, ncol = 1)

# Consolidate legends
legend <- get_legend(legend_2021)

# Arrange all plots with a single legend
theme_set(theme_bw())

grid.arrange(
  richness_2021, evenness_2021, count_2021, 
  richness_2022, evenness_2022, count_2022, 
  richness_2023, evenness_2023, count_2023, 
  ncol = 3,
  top = "Comparison of Richness, Evenness, and Count across Years"
  #,bottom = legend
)




### Prep For Stats ####
Abundance_Data2021 <- full_join(Abundance_ID_aboveground, BurnInfo2021, by = "WS") %>% 
  unite("TreatmentSB",c("Treatment","SB"), sep="_") %>% 
  filter(Date == "2021") %>% 
  filter(Plot == c("2","4"))

Abundance_Data2022 <- full_join(Abundance_ID_aboveground, BurnInfo2021, by = "WS") %>% 
  unite("TreatmentSB",c("Treatment","SB"), sep="_") %>% 
  filter(Date == "2022") 

Abundance_Data2023 <- full_join(Abundance_ID_aboveground, BurnInfo2021, by = "WS") %>% 
  unite("TreatmentSB",c("Treatment","SB"), sep="_") %>% 
  filter(Date == "2023") 

#### 2021 Total Count Stats ####
total_counts2 <- full_join(BurnInfo2021, Abundance_Data2021)
total_counts2 <- total_counts2 %>%
  group_by(TreatmentSB, block, Sample) %>%
  summarise(total_count = sum(Count, na.rm = TRUE)) 

TotalcountModel <- #stores the model output into a named list
  lme(total_count ~ TreatmentSB, #relates richness (or any other dependent variable) to treatment (e.g., ABG and PBG or ABG, PBG0, PBG1, PBG2 for years since burning)
      data = total_counts2, #dataset you are analyzing, this must contain all the data (both treatments and all plots)
      random = ~1|block) #this would be where you'd say north or south unit (which should be a variable in the data frame)
anova.lme(TotalcountModel, type='sequential') #this gives you the ANOVA output from the model, where "sequential" tells it to do a type III anova
emmeans(TotalcountModel, pairwise~TreatmentSB, adjust="tukey") #this gives you contrasts (means and confidence intervals) for each possible pairwise comparison of treatments to know whether they are different or not (overlapping confidence intervals means not different)

# numDF denDF  F-value p-value
# (Intercept)     1    56 136.6163  <.0001
# TreatmentSB     3    56   0.2562  0.8566

#### 2022 Total Count Stats ####
total_counts2 <- full_join(BurnInfo2022, Abundance_Data2022)


total_counts2 <- total_counts2 %>%
  group_by(TreatmentSB, block, Sample) %>%
  summarise(total_count = sum(Count, na.rm = TRUE))

TotalcountModel <- #stores the model output into a named list
  lme(total_count ~ TreatmentSB, #relates richness (or any other dependent variable) to treatment (e.g., ABG and PBG or ABG, PBG0, PBG1, PBG2 for years since burning)
      data = total_counts2, #dataset you are analyzing, this must contain all the data (both treatments and all plots)
      random = ~1|block) #this would be where you'd say north or south unit (which should be a variable in the data frame)
anova.lme(TotalcountModel, type='sequential') #this gives you the ANOVA output from the model, where "sequential" tells it to do a type III anova
emmeans(TotalcountModel, pairwise~TreatmentSB, adjust="tukey") #this gives you contrasts (means and confidence intervals) for each possible pairwise comparison of treatments to know whether they are different or not (overlapping confidence intervals means not different)


# numDF denDF   F-value p-value
# (Intercept)     1    59 17.748082  0.0001
# TreatmentSB     3    59  0.705359  0.5526

#### 2023 Total Count Stats ####
total_counts2 <- full_join(BurnInfo2023, Abundance_Data2023) 

total_counts2 <- total_counts2 %>%
  group_by(TreatmentSB, block, Sample) %>%
  summarise(total_count = sum(Count, na.rm = TRUE))

TotalcountModel <- #stores the model output into a named list
  lme(total_count ~ TreatmentSB, #relates richness (or any other dependent variable) to treatment (e.g., ABG and PBG or ABG, PBG0, PBG1, PBG2 for years since burning)
      data = total_counts2, #dataset you are analyzing, this must contain all the data (both treatments and all plots)
      random = ~1|block) #this would be where you'd say north or south unit (which should be a variable in the data frame)
anova.lme(TotalcountModel, type='sequential') #this gives you the ANOVA output from the model, where "sequential" tells it to do a type III anova
emmeans(TotalcountModel, pairwise~TreatmentSB, adjust="tukey") #this gives you contrasts (means and confidence intervals) for each possible pairwise comparison of treatments to know whether they are different or not (overlapping confidence intervals means not different)

# numDF denDF  F-value p-value
# (Intercept)     1    59 8.965273  0.0040
# TreatmentSB     3    59 2.979528  0.0386

#### 2021 Richness Stats ####

Joined <- full_join(BurnInfo2021, Joined2021) %>% unite("TreatmentSB",c("Treatment","SB"), sep="_")

richModel <- #stores the model output into a named list
  lme(richness ~ TreatmentSB, #relates richness (or any other dependent variable) to treatment (e.g., ABG and PBG or ABG, PBG0, PBG1, PBG2 for years since burnign)
      data = Joined, #dataset you are analyzing, this must contain all the data (both treatments and all plots)
      random = ~1|block) #this would be where you'd say north or south unit (which should be a variable in the dataframe)
anova.lme(richModel, type='sequential') #this gives you the ANOVA output from the model, where "sequential" tells it to do a type III anova
emmeans(richModel, pairwise~TreatmentSB, adjust="tukey") #this gives you contrasts (means and confidence intervals) for each possible pairwise comparison of treatments to know whether they are different or not (overlapping confidence intervals means not different)

# numDF denDF   F-value p-value
# (Intercept)     1    56 197.02313  <.0001
# TreatmentSB     3    56   0.36174  0.7809

#### 2022 Richness Stats ####

Joined <- full_join(BurnInfo2022, Joined2022) %>% unite("TreatmentSB",c("Treatment","SB"), sep="_") %>% na.omit()

richModel <- #stores the model output into a named list
  lme(richness ~ TreatmentSB, #relates richness (or any other dependent variable) to treatment (e.g., ABG and PBG or ABG, PBG0, PBG1, PBG2 for years since burnign)
      data = Joined, #dataset you are analyzing, this must contain all the data (both treatments and all plots)
      random = ~1|block) #this would be where you'd say north or south unit (which should be a variable in the dataframe)
anova.lme(richModel, type='sequential') #this gives you the ANOVA output from the model, where "sequential" tells it to do a type III anova

# numDF denDF  F-value p-value
# (Intercept)     1    59 57.72626  <.0001
# TreatmentSB     3    59  2.57631  0.0623

#### 2023 Richness Stats ####

Joined <- full_join(BurnInfo2023, Joined2023) %>% unite("TreatmentSB",c("Treatment","SB"), sep="_") %>% na.omit()

richModel <- #stores the model output into a named list
  lme(richness ~ TreatmentSB, #relates richness (or any other dependent variable) to treatment (e.g., ABG and PBG or ABG, PBG0, PBG1, PBG2 for years since burnign)
      data = Joined, #dataset you are analyzing, this must contain all the data (both treatments and all plots)
      random = ~1|block) #this would be where you'd say north or south unit (which should be a variable in the dataframe)
anova.lme(richModel, type='sequential') #this gives you the ANOVA output from the model, where "sequential" tells it to do a type III anova
emmeans(richModel, pairwise~TreatmentSB, adjust="tukey") #this gives you contrasts (means and confidence intervals) for each possible pairwise comparison of treatments to know whether they are different or not (overlapping confidence intervals means not different)

# numDF denDF          F-value      p-value
# (Intercept)     1    58 7.863416  0.0069
# TreatmentSB     3    58 4.376781  0.0076

#### 2021 Evenness Stats ####

Joined <- full_join(BurnInfo2021, Joined2021) %>% unite("TreatmentSB",c("Treatment","SB"), sep="_") %>% na.omit()

evarModel <- #stores the model output into a named list
  lme(Evar ~ TreatmentSB, #relates richness (or any other dependent variable) to treatment (e.g., ABG and PBG or ABG, PBG0, PBG1, PBG2 for years since burnign)
      data = Joined, #dataset you are analyzing, this must contain all the data (both treatments and all plots)
      random = ~1|block) #this would be where you'd say north or south unit (which should be a variable in the dataframe)
anova.lme(evarModel, type='sequential') #this gives you the ANOVA output from the model, where "sequential" tells it to do a type III anova

# numDF denDF   F-value p-value
# (Intercept)     1    56 213.94539  <.0001
# TreatmentSB     3    56   2.26334   0.091

#### 2022 Evenness Stats ####

Joined <- full_join(BurnInfo2022, Joined2022) %>% unite("TreatmentSB",c("Treatment","SB"), sep="_") %>% na.omit()

evarModel <- #stores the model output into a named list
  lme(Evar ~ TreatmentSB, #relates richness (or any other dependent variable) to treatment (e.g., ABG and PBG or ABG, PBG0, PBG1, PBG2 for years since burnign)
      data = Joined, #dataset you are analyzing, this must contain all the data (both treatments and all plots)
      random = ~1|block) #this would be where you'd say north or south unit (which should be a variable in the dataframe)
anova.lme(evarModel, type='sequential') #this gives you the ANOVA output from the model, where "sequential" tells it to do a type III anova


# numDF denDF  F-value p-value
# (Intercept)     1    59 2016.432  <.0001
# TreatmentSB     3    59    0.901  0.4462

#### 2023 Evenness Stats ####

Joined <- full_join(BurnInfo2023, Joined2023) %>% unite("TreatmentSB",c("Treatment","SB"), sep="_") %>% na.omit()

evarModel <- #stores the model output into a named list
  lme(Evar ~ TreatmentSB, #relates richness (or any other dependent variable) to treatment (e.g., ABG and PBG or ABG, PBG0, PBG1, PBG2 for years since burnign)
      data = Joined, #dataset you are analyzing, this must contain all the data (both treatments and all plots)
      random = ~1|block) #this would be where you'd say north or south unit (which should be a variable in the dataframe)
anova.lme(evarModel, type='sequential') #this gives you the ANOVA output from the model, where "sequential" tells it to do a type III anova
emmeans(evarModel, pairwise~TreatmentSB, adjust="tukey") #this gives you contrasts (means and confidence intervals) for each possible pairwise comparison of treatments to know whether they are different or not (overlapping confidence intervals means not different)

# numDF denDF   F-value p-value
# (Intercept)     1    58 158.02196  <.0001
# TreatmentSB     3    58   2.46117  0.0716


#### 2021 Weight Stats ####

Joined <- full_join(BurnInfo2021, Z_total_weight_2021) %>% unite("TreatmentSB",c("Treatment","SB"), sep="_") %>% 
  na.omit() %>% 
  mutate(block = ifelse(grepl("S", WS), "North", "South")) %>% 
  filter(Plot %in% c("2", "4"))

weightModel <- #stores the model output into a named list
  lme(combined_weight ~ TreatmentSB, #relates richness (or any other dependent variable) to treatment (e.g., ABG and PBG or ABG, PBG0, PBG1, PBG2 for years since burnign)
      data = Joined, #dataset you are analyzing, this must contain all the data (both treatments and all plots)
      random = ~1|block) #this would be where you'd say north or south unit (which should be a variable in the dataframe)
anova.lme(weightModel, type='sequential') #this gives you the ANOVA output from the model, where "sequential" tells it to do a type III anova

#### 2022 Weight Stats ####

Joined <- full_join(BurnInfo2022, Z_total_weight_2022) %>% unite("TreatmentSB",c("Treatment","SB"), sep="_") %>% 
  na.omit() %>% 
  mutate(block = ifelse(grepl("S", WS), "North", "South"))

weightModel <- #stores the model output into a named list
  lme(combined_weight ~ TreatmentSB, #relates richness (or any other dependent variable) to treatment (e.g., ABG and PBG or ABG, PBG0, PBG1, PBG2 for years since burnign)
      data = Joined, #dataset you are analyzing, this must contain all the data (both treatments and all plots)
      random = ~1|block) #this would be where you'd say north or south unit (which should be a variable in the dataframe)
anova.lme(weightModel, type='sequential') #this gives you the ANOVA output from the model, where "sequential" tells it to do a type III anova

#### 2023 Weight Stats ####

Joined <- full_join(BurnInfo2023, Z_total_weight_2023) %>% unite("TreatmentSB",c("Treatment","SB"), sep="_") %>% 
  na.omit() %>% 
  mutate(block = ifelse(grepl("S", WS), "North", "South"))

weightModel <- #stores the model output into a named list
  lme(combined_weight ~ TreatmentSB, #relates richness (or any other dependent variable) to treatment (e.g., ABG and PBG or ABG, PBG0, PBG1, PBG2 for years since burnign)
      data = Joined, #dataset you are analyzing, this must contain all the data (both treatments and all plots)
      random = ~1|block) #this would be where you'd say north or south unit (which should be a variable in the dataframe)
anova.lme(weightModel, type='sequential') #this gives you the ANOVA output from the model, where "sequential" tells it to do a type III anova


#### CV Richness Calculation 2021 ####
# Getting coefficient of variation (CV) per iteration for richness in 2021

commMetrics2021 <- Abundance_Data2021 %>% 
  filter(TreatmentSB == "ABG_0") %>% 
  community_structure(abundance.var='Count', replicate.var='Sample') 


PBG_average_cv_2021 <- PBG_plot_master_2021 %>%
  group_by(iteration) %>%
  summarize(CV_richness = sd(richness) / mean(richness, na.rm = TRUE))

#Take the mean of the mean for CV PBG richness
PBG_mean_mean_richness_CV_2021 <- mean(PBG_average_cv_2021$CV_richness, na.rm = TRUE)

#ABG Mean
ABG_mean_CV_richness_2021 <- commMetrics2021 %>% 
  summarize(CV_richness = sd(richness) / mean(richness, na.rm = TRUE))


Z_R_CV_2021 <- ((ABG_mean_CV_richness_2021$CV_richness) - (PBG_mean_mean_richness_CV_2021))/(sd(PBG_average_cv_2021$CV_richness))
Z_R_CV_2021
#0.1931119

p_value_R_CV_2021 <- 2*pnorm(-abs(Z_R_CV_2021))

print(p_value_R_CV_2021)
#0.8468713

#### CV Richness Calculation 2022 ####
commMetrics2022 <- Abundance_Data2022 %>%  
  filter(TreatmentSB == "ABG_0") %>% 
  community_structure(abundance.var='Count', replicate.var='Sample')

PBG_average_cv_2022 <- PBG_plot_master_2022 %>%
  group_by(iteration) %>%
  summarize(CV_richness = sd(richness) / mean(richness, na.rm = TRUE))

# Take the mean of the mean for CV PBG richness
PBG_mean_mean_richness_CV_2022 <- mean(PBG_average_cv_2022$CV_richness, na.rm = TRUE)

# ABG Mean
ABG_mean_CV_richness_2022 <- commMetrics2022 %>% 
  summarize(CV_richness = sd(richness) / mean(richness, na.rm = TRUE))

Z_R_CV_2022 <- ((ABG_mean_CV_richness_2022$CV_richness) - (PBG_mean_mean_richness_CV_2022))/(sd(PBG_average_cv_2022$CV_richness))
Z_R_CV_2022
#0.5620573

p_value_R_CV_2022 <- 2*pnorm(-abs(Z_R_CV_2022))

print(p_value_R_CV_2022)
#0.574077

#### CV Richness Calculation 2023 ####
commMetrics2023 <- Abundance_Data2023 %>%
  filter(TreatmentSB == "ABG_0") %>% 
  community_structure(abundance.var='Count', replicate.var='Sample')

PBG_average_cv_2023 <- PBG_plot_master_2023 %>%
  group_by(iteration) %>%
  summarize(CV_richness = sd(richness) / mean(richness, na.rm = TRUE))

# Take the mean of the mean for CV PBG richness
PBG_mean_mean_richness_CV_2023 <- mean(PBG_average_cv_2023$CV_richness, na.rm = TRUE)

# ABG Mean
ABG_mean_CV_richness_2023 <- commMetrics2023 %>% 
  summarize(CV_richness = sd(richness) / mean(richness, na.rm = TRUE))

Z_R_CV_2023 <- ((ABG_mean_CV_richness_2023$CV_richness) - (PBG_mean_mean_richness_CV_2023))/(sd(PBG_average_cv_2023$CV_richness))
Z_R_CV_2023
#0.2951404

p_value_R_CV_2023 <- 2*pnorm(-abs(Z_R_CV_2023))

print(p_value_R_CV_2023)
#0.7678866

#### Richness CV Graph 2021 ####
# Create a data frame for PBG mean CV richness and Z-score
pbg_data <- data.frame(
  metric = "CV_richness",
  CV_richness = PBG_mean_mean_richness_CV_2021,
  treatment = "PBG"
)

# # Create a data frame for ABG mean CV richness
# abg_data <- data.frame(
#   metric = "CV_richness",
#   CV_richness = ABG_mean_CV_richness_2021$CV_richness,
#   treatment = "ABG"
# )

# Combine the data frames
# combined_data <- rbind(pbg_data, abg_data)

# Print the combined data frame
# print(combined_data)


#Richness means graph
CV_Richness_2021 <- ggplot(PBG_average_cv_2021, aes(x = CV_richness, y = ..scaled.., color = "PBG")) +
  geom_density(alpha = 0.5) +
  geom_vline(aes(xintercept = ABG_mean_CV_richness_2021$CV_richness , color = "ABG"), linetype = "dashed", size = 1) +
  labs(title = "2021: Density Plot of Mean CV Richness",
       x = "Mean CV Richness",
       y = "Density") +
  scale_color_manual(values = c("blue", "red"), name = "Legend") +
# annotate("text", x = 20, y = 1, label = "p-value: 0.217", color = "blue", size = 4, hjust = 1, vjust = 1) +
# annotate("text", x = 20, y = 0.80, label = "z-score:  1.24", color = "red", size = 4, hjust = 1, vjust = 1) +
theme(legend.position=c(0.15,0.7),  plot.title = element_text(size = 15),  axis.text.x = element_text(size = 22),  # Adjust the size as needed
      axis.text.y = element_text(size = 22)) 
# scale_x_continuous(limits = c(5, 20)) +
#  scale_y_continuous(limits = c(0, 1)) 

CV_Richness_2021
#CV Count Graph 


#### Richness CV Graph 2022 ####
# Create a data frame for PBG mean CV richness and Z-score
pbg_data <- data.frame(
  metric = "CV_richness",
  CV_richness = PBG_mean_mean_richness_CV_2022,
  treatment = "PBG"
)

# Create a data frame for ABG mean CV richness
# abg_data <- data.frame(
#   metric = "CV_richness",
#   CV_richness = ABG_mean_CV_richness_2022$CV_richness,
#   treatment = "ABG"
# )
# 
# # Combine the data frames
# combined_data <- rbind(pbg_data, abg_data)

# Print the combined data frame
print(combined_data)
#Richness means graph
CV_Richness_2022 <- ggplot(PBG_average_cv_2022, aes(x = CV_richness, y = ..scaled.., color = "PBG")) +
  geom_density(alpha = 0.5) +
  geom_vline(aes(xintercept = ABG_mean_CV_richness_2022$CV_richness , color = "ABG"), linetype = "dashed", size = 1) +
  labs(title = "2022: Density Plot of Mean CV Richness",
       x = "Mean CV Richness",
       y = "Density") +
  scale_color_manual(values = c("blue", "red"), name = "Legend") +
  # annotate("text", x = 20, y = 1, label = "p-value: 0.217", color = "blue", size = 4, hjust = 1, vjust = 1) +
  # annotate("text", x = 20, y = 0.80, label = "z-score:  1.24", color = "red", size = 4, hjust = 1, vjust = 1) +
  theme(legend.position = "none",  plot.title = element_text(size = 15),  axis.text.x = element_text(size = 22),  # Adjust the size as needed
        axis.text.y = element_text(size = 22)) 
# scale_x_continuous(limits = c(5, 20)) +
#  scale_y_continuous(limits = c(0, 1)) 

CV_Richness_2022

#CV Count Graph 
#### Richness CV Graph 2023 ####
# Create a data frame for PBG mean CV richness and Z-score
pbg_data <- data.frame(
  metric = "CV_richness",
  CV_richness = PBG_mean_mean_richness_CV_2023,
  treatment = "PBG"
)

# # Create a data frame for ABG mean CV richness
# abg_data <- data.frame(
#   metric = "CV_richness",
#   CV_richness = ABG_mean_CV_richness_2023$CV_richness,
#   treatment = "ABG"
# )
# 
# # Combine the data frames
# combined_data <- rbind(pbg_data, abg_data)

# Print the combined data frame
print(combined_data)

# CV Count Graph 
CV_Richness_2023 <- ggplot(PBG_average_cv_2023, aes(x = CV_richness, y = ..scaled.., color = "PBG")) +
  geom_density(alpha = 0.5) +
  geom_vline(aes(xintercept = ABG_mean_CV_richness_2023$CV_richness , color = "ABG"), linetype = "dashed", size = 1) +
  labs(title = "2023: Density Plot of Mean CV Richness",
       x = "Mean CV Richness",
       y = "Density") +
  scale_color_manual(values = c("blue", "red"), name = "Legend") +
  theme(legend.position = "none",  plot.title = element_text(size = 15),  axis.text.x = element_text(size = 22),  # Adjust the size as needed
        axis.text.y = element_text(size = 22)) 
# annotate("text", x = 20, y = 1, label = "p-value: 0.217", color = "blue", size = 4, hjust = 1, vjust = 1) +
# annotate("text", x = 20, y = 0.80, label = "z-score:  1.24", color = "red", size = 4, hjust = 1, vjust = 1) +
# scale_x_continuous(limits = c(5, 20)) +
# scale_y_continuous(limits = c(0, 1)) 

CV_Richness_2023

#### CV Evar Calculation 2021 ####
commMetrics2021 <- Abundance_Data2021 %>%  
  filter(TreatmentSB == "ABG_0") %>% 
  community_structure(abundance.var='Count', replicate.var='Sample')

PBG_average_cv_2021 <- PBG_plot_master_2021 %>%
  group_by(iteration) %>%
  summarize(CV_Evar = sd(Evar) / mean(Evar, na.rm = TRUE))

# Take the mean of the mean for CV PBG Evar
PBG_mean_mean_Evar_CV_2021 <- mean(PBG_average_cv_2021$CV_Evar, na.rm = TRUE)

# ABG Mean
ABG_mean_CV_Evar_2021 <- commMetrics2021 %>% 
  summarize(CV_Evar = sd(Evar) / mean(Evar, na.rm = TRUE))

Z_E_CV_2021 <- ((ABG_mean_CV_Evar_2021$CV_Evar) - (PBG_mean_mean_Evar_CV_2021))/(sd(PBG_average_cv_2021$CV_Evar))
Z_E_CV_2021
#2.019622

p_value_E_CV_2021 <- 2*pnorm(-abs(Z_E_CV_2021))

print(p_value_E_CV_2021)
#0.04342265

#### CV Evar Calculation 2022 ####
commMetrics2022 <- Abundance_Data2022 %>%  
  filter(TreatmentSB == "ABG_0") %>% 
  community_structure(abundance.var='Count', replicate.var='Sample')

PBG_average_cv_2022 <- PBG_plot_master_2022 %>%
  group_by(iteration) %>%
  summarize(CV_Evar = sd(Evar) / mean(Evar, na.rm = TRUE))

# Take the mean of the mean for CV PBG Evar
PBG_mean_mean_Evar_CV_2022 <- mean(PBG_average_cv_2022$CV_Evar, na.rm = TRUE)

# ABG Mean
ABG_mean_CV_Evar_2022 <- commMetrics2022 %>% 
  summarize(CV_Evar = sd(Evar, na.rm = T) / mean(Evar, na.rm = TRUE))

Z_E_CV_2022 <- ((ABG_mean_CV_Evar_2022$CV_Evar) - (PBG_mean_mean_Evar_CV_2022))/(sd(PBG_average_cv_2022$CV_Evar))
Z_E_CV_2022
#-0.2505794

p_value_E_CV_2022 <- 2*pnorm(-abs(Z_E_CV_2022))

print(p_value_E_CV_2022)
#0.8021393


#### CV Evar Calculation 2023 ####

commMetrics2023 <- Abundance_Data2023 %>% 
  filter(TreatmentSB == "ABG_0") %>% 
  community_structure(abundance.var='Count', replicate.var='Sample')

PBG_average_cv_2023 <- PBG_plot_master_2023 %>%
  group_by(iteration) %>%
  summarize(CV_Evar = sd(Evar) / mean(Evar, na.rm = TRUE))

# Take the mean of the mean for CV PBG Evar
PBG_mean_mean_Evar_CV_2023 <- mean(PBG_average_cv_2023$CV_Evar, na.rm = TRUE)

# ABG Mean
ABG_mean_CV_Evar_2023 <- commMetrics2023 %>% 
  summarize(CV_Evar = sd(Evar, na.rm = TRUE) / mean(Evar, na.rm = TRUE))

Z_E_CV_2023 <- ((ABG_mean_CV_Evar_2023$CV_Evar) - (PBG_mean_mean_Evar_CV_2023))/(sd(PBG_average_cv_2023$CV_Evar))
Z_E_CV_2023
#1.029132

p_value_E_CV_2023 <- 2*pnorm(-abs(Z_E_CV_2023))

print(p_value_E_CV_2023)
#0.3034176

#### Evar CV Graph 2021 ####
# Create a data frame for PBG mean CV Evar and Z-score
pbg_data <- data.frame(
  metric = "CV_Evar",
  CV_Evar = PBG_mean_mean_Evar_CV_2021,
  treatment = "PBG"
)

# Create a data frame for ABG mean CV Evar
# abg_data <- data.frame(
#   metric = "CV_Evar",
#   CV_Evar = ABG_mean_CV_Evar_2021$CV_Evar,
#   treatment = "ABG"
# )

# Combine the data frames
# combined_data <- rbind(pbg_data, abg_data)

# Print the combined data frame
print(combined_data)

# CV Evar Graph 
CV_Evar_2021 <- ggplot(PBG_average_cv_2021, aes(x = CV_Evar, y = ..scaled.., color = "PBG")) +
  geom_density(alpha = 0.5) +
  geom_vline(aes(xintercept = ABG_mean_CV_Evar_2021$CV_Evar , color = "ABG"), linetype = "dashed", size = 1) +
  labs(title = "2021: Density Plot of Mean CV Evenness",
       x = "Mean CV Evenness",
       y = "Density") +
  scale_color_manual(values = c("blue", "red"), name = "Legend") +
# annotate("text", x = 20, y = 1, label = "p-value: 0.217", color = "blue", size = 4, hjust = 1, vjust = 1) +
# annotate("text", x = 20, y = 0.80, label = "z-score:  1.24", color = "red", size = 4, hjust = 1, vjust = 1) +
theme(legend.position = "none",  plot.title = element_text(size = 15),  axis.text.x = element_text(size = 22),  # Adjust the size as needed
      axis.text.y = element_text(size = 22)) 
# scale_x_continuous(limits = c(5, 20)) +
# scale_y_continuous(limits = c(0, 1)) 

CV_Evar_2021

#### Evar CV Graph 2022 ####
# Create a data frame for PBG mean CV Evar and Z-score
pbg_data <- data.frame(
  metric = "CV_Evar",
  CV_Evar = PBG_mean_mean_Evar_CV_2022,
  treatment = "PBG"
)

# Create a data frame for ABG mean CV Evar
# abg_data <- data.frame(
#   metric = "CV_Evar",
#   CV_Evar = ABG_mean_CV_Evar_2022$CV_Evar,
#   treatment = "ABG"
# )
# 
# # Combine the data frames
# combined_data <- rbind(pbg_data, abg_data)

# Print the combined data frame
# print(combined_data)

# CV Evar Graph 
CV_Evar_2022 <- ggplot(PBG_average_cv_2022, aes(x = CV_Evar, y = ..scaled.., color = "PBG")) +
  geom_density(alpha = 0.5) +
  geom_vline(aes(xintercept = ABG_mean_CV_Evar_2022$CV_Evar , color = "ABG"), linetype = "dashed", size = 1) +
  labs(title = "2022: Density Plot of Mean CV Evenness",
       x = "Mean CV Evenness",
       y = "Density") +
  scale_color_manual(values = c("blue", "red"), name = "Legend") +
# annotate("text", x = 20, y = 1, label = "p-value: 0.217", color = "blue", size = 4, hjust = 1, vjust = 1) +
# annotate("text", x = 20, y = 0.80, label = "z-score:  1.24", color = "red", size = 4, hjust = 1, vjust = 1) +
 theme(legend.position = "none",  plot.title = element_text(size = 15),  axis.text.x = element_text(size = 22),  # Adjust the size as needed
       axis.text.y = element_text(size = 22)) 
# scale_x_continuous(limits = c(5, 20)) +
# scale_y_continuous(limits = c(0, 1)) 

CV_Evar_2022


#### Evar CV Graph 2023 ####
# Create a data frame for PBG mean CV Evar and Z-score
pbg_data <- data.frame(
  metric = "CV_Evar",
  CV_Evar = PBG_mean_mean_Evar_CV_2023,
  treatment = "PBG"
)

# Create a data frame for ABG mean CV Evar
# abg_data <- data.frame(
#   metric = "CV_Evar",
#   CV_Evar = ABG_mean_CV_Evar_2023$CV_Evar,
#   treatment = "ABG"
# )
# 
# # Combine the data frames
# combined_data <- rbind(pbg_data, abg_data)
# 
# # Print the combined data frame
# print(combined_data)

# CV Evar Graph 
CV_Evar_2023 <- ggplot(PBG_average_cv_2023, aes(x = CV_Evar, y = ..scaled.., color = "PBG")) +
  geom_density(alpha = 0.5) +
  geom_vline(aes(xintercept = ABG_mean_CV_Evar_2023$CV_Evar , color = "ABG"), linetype = "dashed", size = 1) +
  labs(title = "2023: Density Plot of Mean CV Evenness",
       x = "Mean CV Evenness",
       y = "Density") +
  scale_color_manual(values = c("blue", "red"), name = "Legend") +
# annotate("text", x = 20, y = 1, label = "p-value: 0.217", color = "blue", size = 4, hjust = 1, vjust = 1) +
# annotate("text", x = 20, y = 0.80, label = "z-score:  1.24", color = "red", size = 4, hjust = 1, vjust = 1) +
theme(legend.position = "none",  plot.title = element_text(size = 15),  axis.text.x = element_text(size = 22),  # Adjust the size as needed
      axis.text.y = element_text(size = 22))
# scale_x_continuous(limits = c(5, 20)) +
# scale_y_continuous(limits = c(0, 1)) 

CV_Evar_2023



#### CV Count Calculation 2021 ####
countMetrics2021 <- Abundance_Data2021



PBG_average_cv_2021 <- PBG_plot_master_2021 %>%
  group_by(iteration) %>%
  summarize(CV_count = sd(Count, na.rm = T) / mean(Count, na.rm = TRUE))

#Take the mean of the mean for CV PBG count
PBG_mean_mean_count_CV_2021 <- mean(PBG_average_cv_2021$CV_count, na.rm = TRUE)

#ABG Mean
ABG_mean_CV_count_2021 <- countMetrics2021 %>%
  filter(TreatmentSB == "ABG_0") %>%
  summarize(CV_count = sd(Count) / mean(Count, na.rm = TRUE))


Z_C_CV_2021 <- ((ABG_mean_CV_count_2021$CV_count) - (PBG_mean_mean_count_CV_2021))/(sd(PBG_average_cv_2021$CV_count))
Z_C_CV_2021
#10.92508

p_value_C_CV_2021 <- 2*pnorm(-abs(Z_C_CV_2021))

print(p_value_C_CV_2021)
# 8.746414e-28
#### CV Count Calculation 2022 ####
countMetrics2022 <- Abundance_Data2022 %>% 
  filter(TreatmentSB == "ABG_0")

PBG_average_cv_2022 <- PBG_plot_master_2022 %>%
  group_by(iteration) %>%
  summarize(CV_count = sd(Count, na.rm = TRUE) / mean(Count, na.rm = TRUE))

# Take the mean of the mean for CV PBG count
PBG_mean_mean_count_CV_2022 <- mean(PBG_average_cv_2022$CV_count, na.rm = TRUE)

# ABG Mean
ABG_mean_CV_count_2022 <- countMetrics2022 %>% 
  summarize(CV_count = sd(Count) / mean(Count, na.rm = TRUE))

Z_C_CV_2022 <- ((ABG_mean_CV_count_2022$CV_count) - (PBG_mean_mean_count_CV_2022))/(sd(PBG_average_cv_2022$CV_count))
Z_C_CV_2022
#7.042279

p_value_C_CV_2022 <- 2*pnorm(-abs(Z_C_CV_2022))

print(p_value_C_CV_2022)
#1.891211e-12

#### CV Count Calculation 2023 ####
countMetrics2023 <- Abundance_Data2023 %>% 
  filter(TreatmentSB == "ABG_0")

PBG_average_cv_2023 <- PBG_plot_master_2023 %>%
  group_by(iteration) %>%
  summarize(CV_count = sd(Count, na.rm = TRUE) / mean(Count, na.rm = TRUE))

# Take the mean of the mean for CV PBG count
PBG_mean_mean_count_CV_2023 <- mean(PBG_average_cv_2023$CV_count, na.rm = TRUE)

# ABG Mean
ABG_mean_CV_count_2023 <- countMetrics2023 %>% 
  summarize(CV_count = sd(Count) / mean(Count, na.rm = TRUE))

Z_C_CV_2023 <- ((ABG_mean_CV_count_2023$CV_count) - (PBG_mean_mean_count_CV_2023))/(sd(PBG_average_cv_2023$CV_count))
Z_C_CV_2023
#6.090173

p_value_C_CV_2023 <- 2*pnorm(-abs(Z_C_CV_2023))

print(p_value_C_CV_2023)
#1.12789e-09
#### Count CV Graph 2021 ####
# Create a data frame for PBG mean CV Count and Z-score
pbg_data <- data.frame(
  metric = "CV_Count",
  CV_Count = PBG_mean_mean_count_CV_2021,
  treatment = "PBG"
)

# Create a data frame for ABG mean CV Count
abg_data <- data.frame(
  metric = "CV_Count",
  CV_Count = ABG_mean_CV_count_2021$CV_count,
  treatment = "ABG"
)

# Combine the data frames
combined_data <- rbind(pbg_data, abg_data)

# Print the combined data frame
print(combined_data)

# CV Count Graph 
CV_Count_2021 <- ggplot(PBG_average_cv_2021, aes(x = CV_count, y = ..scaled.., color = "PBG")) +
  geom_density(alpha = 0.5) +
  geom_vline(aes(xintercept = ABG_mean_CV_count_2021 $CV_count , color = "ABG"), linetype = "dashed", size = 1) +
  labs(title = "2021: Density Plot of Mean CV Abundance",
       x = "Mean CV Abundance",
       y = "Density") +
  scale_color_manual(values = c("blue", "red"), name = "Legend") +
# annotate("text", x = 20, y = 1, label = "p-value: 0.217", color = "blue", size = 4, hjust = 1, vjust = 1) +
# annotate("text", x = 20, y = 0.80, label = "z-score:  1.24", color = "red", size = 4, hjust = 1, vjust = 1) +
theme(legend.position = "none",  plot.title = element_text(size = 15),  axis.text.x = element_text(size = 22),  # Adjust the size as needed
      axis.text.y = element_text(size = 22)) + 
scale_x_continuous(limits = c(0.2, 1.6)) +
scale_y_continuous(limits = c(0, 1)) 

CV_Count_2021


#### Count CV Graph 2022 ####
# Create a data frame for PBG mean CV Count and Z-score
pbg_data <- data.frame(
  metric = "CV_Count",
  CV_Count = PBG_mean_mean_count_CV_2022,
  treatment = "PBG"
)

# Create a data frame for ABG mean CV Count
abg_data <- data.frame(
  metric = "CV_Count",
  CV_Count = ABG_mean_CV_count_2022$CV_count,
  treatment = "ABG"
)

# Combine the data frames
combined_data <- rbind(pbg_data, abg_data)

# Print the combined data frame
print(combined_data)

# CV Count Graph 
CV_Count_2022 <- ggplot(PBG_average_cv_2022, aes(x = CV_count, y = ..scaled.., color = "PBG")) +
  geom_density(alpha = 0.5) +
  geom_vline(aes(xintercept = ABG_mean_CV_count_2022$CV_count , color = "ABG"), linetype = "dashed", size = 1) +
  labs(title = "2022: Density Plot of Mean CV Abundance",
       x = "Mean CV Abundance",
       y = "Density") +
  scale_color_manual(values = c("blue", "red"), name = "Legend") +
# annotate("text", x = 20, y = 1, label = "p-value: 0.217", color = "blue", size = 4, hjust = 1, vjust = 1) +
# annotate("text", x = 20, y = 0.80, label = "z-score:  1.24", color = "red", size = 4, hjust = 1, vjust = 1) +
theme(legend.position = "none",  plot.title = element_text(size = 15),  axis.text.x = element_text(size = 22),  # Adjust the size as needed
      axis.text.y = element_text(size = 22)) +
scale_x_continuous(limits = c(0.2, 1.6)) +
scale_y_continuous(limits = c(0, 1)) 

CV_Count_2022


#### Count CV Graph 2023 ####
# Create a data frame for PBG mean CV Count and Z-score
pbg_data <- data.frame(
  metric = "CV_Count",
  CV_Count = PBG_mean_mean_count_CV_2023,
  treatment = "PBG"
)

# Create a data frame for ABG mean CV Count
abg_data <- data.frame(
  metric = "CV_Count",
  CV_Count = ABG_mean_CV_count_2023$CV_count,
  treatment = "ABG"
)

# Combine the data frames
combined_data <- rbind(pbg_data, abg_data)

# Print the combined data frame
print(combined_data)

# CV Count Graph 
CV_Count_2023 <- ggplot(PBG_average_cv_2023, aes(x = CV_count, y = ..scaled.., color = "PBG")) +
  geom_density(alpha = 0.5) +
  geom_vline(aes(xintercept = ABG_mean_CV_count_2023$CV_count , color = "ABG"), linetype = "dashed", size = 1) +
  labs(title = "2023: Density Plot of Mean CV Abundance",
       x = "Mean CV Abundance",
       y = "Density") +
  scale_color_manual(values = c("blue", "red"), name = "Legend") +
# annotate("text", x = 20, y = 1, label = "p-value: 0.217", color = "blue", size = 4, hjust = 1, vjust = 1) +
# annotate("text", x = 20, y = 0.80, label = "z-score:  1.24", color = "red", size = 4, hjust = 1, vjust = 1) +
theme(legend.position = "none",  plot.title = element_text(size = 15),  axis.text.x = element_text(size = 22),  # Adjust the size as needed
      axis.text.y = element_text(size = 22)) +
scale_x_continuous(limits = c(0.2, 1.6)) +
scale_y_continuous(limits = c(0, 1)) 

CV_Count_2023


#### Big CV Graphs ####
multi_panel_graph <- grid.arrange(
  CV_Richness_2021, CV_Evar_2021, CV_Count_2021,
  CV_Richness_2022, CV_Evar_2022, CV_Count_2022,
  CV_Richness_2023, CV_Evar_2023, CV_Count_2023,
  nrow = 3, ncol = 3
)



#### Years Since Burned Graph Settings ####
Axis_Label_Size_1 <- 20 #The actualy ticks and treatments
Axis_Text_Size <- 20 #Text on axislike abuandance, richness etc
#### 2021 Years Since Burned, Count Evenness Graphs ####
# Box plot for count
count_boxplot <- ggplot(Abundance_Data2021, aes(x=TreatmentSB, y=Count, fill=TreatmentSB)) +
  geom_boxplot() +
  labs(title="Abundance Comparison for Treatments with Years Since Burned (2021)", x="Treatment", y="Abundance") +
  scale_fill_manual(values = c("ABG" = "blue", "PBG" = "red", "PBG_0" = "red", "PBG_1" = "red", "PBG_2" = "red", "ABG_0" = "blue")) + # Color PBG_0, PBG_1, PBG_2 as red, ABG_0 as blue
  theme_minimal() +
  #annotate("text", x = -Inf, y = Inf, label = "P value: 0.857", vjust = 1, hjust = 0, size = 3.5, color = "black") + 
  theme(legend.position = "none", plot.title = element_text(size = 9),
        axis.text.x = element_text(size = Axis_Label_Size_1), # Increase the size of x-axis labels
        axis.text.y = element_text(size = Axis_Label_Size_1),
        axis.title.x = element_text(size = Axis_Label_Size_1), # Increase the size of x-axis title
        axis.title.y = element_text(size = Axis_Label_Size_1)) + # Increase the size of y-axis title
  scale_x_discrete(labels = c("ABG_0" = "ABG 0", "PBG_0" = "PBG 0", "PBG_1" = "PBG 1", "PBG_2" = "PBG 2")) + # Set custom axis labels
  theme(legend.position = "none", plot.title = element_text(size = 9),
        axis.text.x = element_text(angle = 0, hjust = 0.5)) # Adjust x-axis text angle for

# Display the plot
count_boxplot


# numDF denDF  F-value p-value
# (Intercept)     1    56 136.6163  <.0001
# TreatmentSB     3    56   0.2562  0.8566

# Boxplot for Richness
richness_boxplot <- ggplot(Joined2021, aes(x=TreatmentSB, y=richness, fill=TreatmentSB)) +
  geom_boxplot() +
  labs(title="Richness Comparison for Treatments with Years Since Burned (2021)", x="Treatment", y="Richness") +
  scale_fill_manual(values = c("ABG_0" = "blue", "PBG_0" = "red", "PBG_1" = "red", "PBG_2" = "red")) + # Color PBG_0, PBG_1, PBG_2 as red, ABG_0 as blue
  theme_minimal() +
 # annotate("text", x = -Inf, y = Inf, label = "P value: 0.725", vjust = 1, hjust = 0, size = 3.5, color = "black") + 
  theme(legend.position = "none", plot.title = element_text(size = 9),
        axis.text.x = element_text(size = Axis_Label_Size_1), # Increase the size of x-axis labels
        axis.text.y = element_text(size = Axis_Label_Size_1),
        axis.title.x = element_text(size = Axis_Label_Size_1), # Increase the size of x-axis title
        axis.title.y = element_text(size = Axis_Label_Size_1)) + # Increase the size of y-axis title
  scale_x_discrete(labels = c("ABG_0" = "ABG 0", "PBG_0" = "PBG 0", "PBG_1" = "PBG 1", "PBG_2" = "PBG 2")) + # Set custom axis labels
  theme(legend.position = "none", plot.title = element_text(size = 9),
        axis.text.x = element_text(angle = 0, hjust = 0.5)) # Adjust x-axis text angle for
# Display the plot
richness_boxplot

# numDF denDF  F-value p-value
# (Intercept)     1    56 52.75320  <.0001
# TreatmentSB     3    56  0.43999  0.7253


# Boxplot for Evenness
evenness_boxplot <- ggplot(Joined2021, aes(x=TreatmentSB, y=Evar, fill=TreatmentSB)) +
  geom_boxplot() +
  labs(title="Evenness Comparison for Treatments with Years Since Burned (2021)", x="Treatment", y="Evenness") +
  scale_fill_manual(values = c("ABG_0" = "blue", "PBG_0" = "red", "PBG_1" = "red", "PBG_2" = "red")) + # Color PBG_0, PBG_1, PBG_2 as red, ABG_0 as blue
  theme_minimal() +
  #annotate("text", x = -Inf, y = Inf, label = "P value: 0.091", vjust = 1, hjust = 0, size = 3.5, color = "black") + 
  theme(legend.position = "none", plot.title = element_text(size = 9),
        axis.text.x = element_text(size = Axis_Label_Size_1), # Increase the size of x-axis labels
        axis.text.y = element_text(size = Axis_Label_Size_1),
        axis.title.x = element_text(size = Axis_Label_Size_1), # Increase the size of x-axis title
        axis.title.y = element_text(size = Axis_Label_Size_1)) + # Increase the size of y-axis title
  scale_x_discrete(labels = c("ABG_0" = "ABG 0", "PBG_0" = "PBG 0", "PBG_1" = "PBG 1", "PBG_2" = "PBG 2")) + # Set custom axis labels
  theme(legend.position = "none", plot.title = element_text(size = 9),
        axis.text.x = element_text(angle = 0, hjust = 0.5)) # Adjust x-axis text angle for
# Display the plot
evenness_boxplot

# numDF denDF   F-value p-value
# (Intercept)     1    56 213.94539  <.0001
# TreatmentSB     3    56   2.26334   0.091


#### 2022 Years Since Burned, Count Evenness Graphs ####
# Box plot for count
count_boxplot_2022 <- ggplot(Abundance_Data2022, aes(x=TreatmentSB, y=Count, fill=TreatmentSB)) +
  geom_boxplot() +
  labs(title="Abundance Comparison for Treatments with Years Since Burned (2022)", x="Treatment", y="Abundance") +
  scale_fill_manual(values = c("ABG" = "blue", "PBG" = "red", "PBG_0" = "red", "PBG_1" = "red", "PBG_2" = "red", "ABG_0" = "blue")) + # Color PBG_0, PBG_1, PBG_2 as red, ABG_0 as blue
  theme_minimal() +
 # annotate("text", x = -Inf, y = Inf, label = "P value: 0.553", vjust = 1, hjust = 0, size = 3.5, color = "black") + 
  theme(legend.position = "none", plot.title = element_text(size = 9),
        axis.text.x = element_text(size = Axis_Label_Size_1), # Increase the size of x-axis labels
        axis.text.y = element_text(size = Axis_Label_Size_1),
        axis.title.x = element_text(size = Axis_Label_Size_1), # Increase the size of x-axis title
        axis.title.y = element_text(size = Axis_Label_Size_1)) + # Increase the size of y-axis title
  scale_x_discrete(labels = c("ABG_0" = "ABG 0", "PBG_0" = "PBG 0", "PBG_1" = "PBG 1", "PBG_2" = "PBG 2")) + # Set custom axis labels
  theme(legend.position = "none", plot.title = element_text(size = 9),
        axis.text.x = element_text(angle = 0, hjust = 0.5)) # Adjust x-axis text angle for

# Display the plot
count_boxplot_2022

# numDF denDF   F-value p-value
# (Intercept)     1    59 17.748082  0.0001
# TreatmentSB     3    59  0.705359  0.5526

# Boxplot for Richness
richness_boxplot_2022 <- ggplot(Joined2022, aes(x=TreatmentSB, y=richness, fill=TreatmentSB)) +
  geom_boxplot() +
  labs(title="Richness Comparison for Treatments with Years Since Burned (2022)", x="Treatment", y="Richness") +
  scale_fill_manual(values = c("ABG_0" = "blue", "PBG_0" = "red", "PBG_1" = "red", "PBG_2" = "red")) + # Color PBG_0, PBG_1, PBG_2 as red, ABG_0 as blue
  theme_minimal() +
 # annotate("text", x = -Inf, y = Inf, label = "P value: 0.062", vjust = 1, hjust = 0, size = 3.5, color = "black") + 
  theme(legend.position = "none", plot.title = element_text(size = 9),
        axis.text.x = element_text(size = Axis_Label_Size_1), # Increase the size of x-axis labels
        axis.text.y = element_text(size = Axis_Label_Size_1),
        axis.title.x = element_text(size = Axis_Label_Size_1), # Increase the size of x-axis title
        axis.title.y = element_text(size = Axis_Label_Size_1)) + # Increase the size of y-axis title
  scale_x_discrete(labels = c("ABG_0" = "ABG 0", "PBG_0" = "PBG 0", "PBG_1" = "PBG 1", "PBG_2" = "PBG 2")) + # Set custom axis labels
  theme(legend.position = "none", plot.title = element_text(size = 9),
        axis.text.x = element_text(angle = 0, hjust = 0.5)) # Adjust x-axis text angle for

# numDF denDF  F-value p-value
# (Intercept)     1    59 57.72626  <.0001
# TreatmentSB     3    59  2.57631  0.0623

# Display the plot
richness_boxplot_2022

# Boxplot for Evenness
evenness_boxplot_2022 <- ggplot(Joined2022, aes(x=TreatmentSB, y=Evar, fill=TreatmentSB)) +
  geom_boxplot() +
  labs(title="Evenness Comparison for Treatments with Years Since Burned (2022)", x="Treatment", y="Evenness") +
  scale_fill_manual(values = c("ABG_0" = "blue", "PBG_0" = "red", "PBG_1" = "red", "PBG_2" = "red")) + # Color PBG_0, PBG_1, PBG_2 as red, ABG_0 as blue
  theme_minimal() +
#  annotate("text", x = -Inf, y = Inf, label = "P value: 0.446", vjust = 1, hjust = 0, size = 3.5, color = "black") + 
  theme(legend.position = "none", plot.title = element_text(size = 9),
        axis.text.x = element_text(size = Axis_Label_Size_1), # Increase the size of x-axis labels
        axis.text.y = element_text(size = Axis_Label_Size_1),
        axis.title.x = element_text(size = Axis_Label_Size_1), # Increase the size of x-axis title
        axis.title.y = element_text(size = Axis_Label_Size_1)) + # Increase the size of y-axis title
  scale_x_discrete(labels = c("ABG_0" = "ABG 0", "PBG_0" = "PBG 0", "PBG_1" = "PBG 1", "PBG_2" = "PBG 2")) + # Set custom axis labels
  theme(legend.position = "none", plot.title = element_text(size = 9),
        axis.text.x = element_text(angle = 0, hjust = 0.5)) # Adjust x-axis text angle for
# Display the plot
evenness_boxplot_2022

# numDF denDF  F-value p-value
# (Intercept)     1    59 2016.432  <.0001
# TreatmentSB     3    59    0.901  0.4462


#### 2023 Years Since Burned, Count Evenness Graphs ####
# Box plot for count
count_boxplot_2023 <- ggplot(Abundance_Data2023, aes(x=TreatmentSB, y=Count, fill=TreatmentSB)) +
  geom_boxplot() +
  labs(title="Abundance Comparison for Treatments with Years Since Burned (2023)", x="Treatment", y="Abundance") +
  scale_fill_manual(values = c("ABG" = "blue", "PBG" = "red", "PBG_0" = "red", "PBG_1" = "red", "PBG_2" = "red", "ABG_0" = "blue")) + # Color PBG_0, PBG_1, PBG_2 as red, ABG_0 as blue
  scale_x_discrete(labels = c("ABG_0" = "ABG 0", "PBG_0" = "PBG 0", "PBG_1" = "PBG 1", "PBG_2" = "PBG 2")) + # Set custom axis labels
  theme_minimal() +
#  annotate("text", x = -Inf, y = Inf, label = "P value: 0.039", vjust = 1, hjust = 0, size = 3.5, color = "black") + 
  theme(legend.position = "none", plot.title = element_text(size = 9),
        axis.text.x = element_text(size = Axis_Label_Size_1), # Increase the size of x-axis labels
        axis.text.y = element_text(size = Axis_Label_Size_1),
        axis.title.x = element_text(size = Axis_Label_Size_1), # Increase the size of x-axis title
        axis.title.y = element_text(size = Axis_Label_Size_1)) + # Increase the size of y-axis title
  scale_x_discrete(labels = c("ABG_0" = "ABG 0", "PBG_0" = "PBG 0", "PBG_1" = "PBG 1", "PBG_2" = "PBG 2")) + # Set custom axis labels
  theme(legend.position = "none", plot.title = element_text(size = 9),
        axis.text.x = element_text(angle = 0, hjust = 0.5)) # Adjust x-axis text angle for
# Display the plot
count_boxplot_2023

# numDF denDF  F-value p-value
# (Intercept)     1    59 8.965273  0.0040
# TreatmentSB     3    59 2.979528  0.0386


# Boxplot for Richness
richness_boxplot_2023 <- ggplot(Joined2023, aes(x=TreatmentSB, y=richness, fill=TreatmentSB)) +
  geom_boxplot() +
  labs(title="Richness Comparison for Treatments with Years Since Burned (2023)", x="Treatment", y="Richness") +
  scale_fill_manual(values = c("ABG_0" = "blue", "PBG_0" = "red", "PBG_1" = "red", "PBG_2" = "red")) + # Color PBG_0, PBG_1, PBG_2 as red, ABG_0 as blue
  scale_x_discrete(labels = c("ABG_0" = "ABG 0", "PBG_0" = "PBG 0", "PBG_1" = "PBG 1", "PBG_2" = "PBG 2")) + # Set custom axis labels
  theme_minimal() +
#  annotate("text", x = -Inf, y = Inf, label = "P value: 0.008", vjust = 1, hjust = 0, size = 3.5, color = "black") + 
  theme(legend.position = "none", plot.title = element_text(size = 9),
        axis.text.x = element_text(size = Axis_Label_Size_1), # Increase the size of x-axis labels
        axis.text.y = element_text(size = Axis_Label_Size_1),
        axis.title.x = element_text(size = Axis_Label_Size_1), # Increase the size of x-axis title
        axis.title.y = element_text(size = Axis_Label_Size_1)) + # Increase the size of y-axis title
  scale_x_discrete(labels = c("ABG_0" = "ABG 0", "PBG_0" = "PBG 0", "PBG_1" = "PBG 1", "PBG_2" = "PBG 2")) + # Set custom axis labels
  theme(legend.position = "none", plot.title = element_text(size = 9),
        axis.text.x = element_text(angle = 0, hjust = 0.5)) # Adjust x-axis text angle for
richness_boxplot_2023

# numDF denDF          F-value      p-value
# (Intercept)     1    58 7.863416  0.0069
# TreatmentSB     3    58 4.376781  0.0076

# Boxplot for Evenness
evenness_boxplot_2023 <- ggplot(Joined2023, aes(x=TreatmentSB, y=Evar, fill=TreatmentSB)) +
  geom_boxplot() +
  labs(title="Evenness Comparison for Treatments with Years Since Burned (2023)", x="Treatment", y="Evenness") +
  scale_fill_manual(values = c("ABG_0" = "blue", "PBG_0" = "red", "PBG_1" = "red", "PBG_2" = "red")) + # Color PBG_0, PBG_1, PBG_2 as red, ABG_0 as blue
  scale_x_discrete(labels = c("ABG_0" = "ABG 0", "PBG_0" = "PBG 0", "PBG_1" = "PBG 1", "PBG_2" = "PBG 2")) + # Set custom axis labels
  theme_minimal() +
# annotate("text", x = -Inf, y = Inf, label = "P value: 0.072", vjust = 1, hjust = 0, size = 3.5, color = "black") + 
  theme(legend.position = "none", plot.title = element_text(size = 9),
        axis.text.x = element_text(size = Axis_Label_Size_1), # Increase the size of x-axis labels
        axis.text.y = element_text(size = Axis_Label_Size_1),
        axis.title.x = element_text(size = Axis_Label_Size_1), # Increase the size of x-axis title
        axis.title.y = element_text(size = Axis_Label_Size_1)) + # Increase the size of y-axis title
  scale_x_discrete(labels = c("ABG_0" = "ABG 0", "PBG_0" = "PBG 0", "PBG_1" = "PBG 1", "PBG_2" = "PBG 2")) + # Set custom axis labels
  theme(legend.position = "none", plot.title = element_text(size = 9),
        axis.text.x = element_text(angle = 0, hjust = 0.5)) # Adjust x-axis text angle for
# Display the plot
evenness_boxplot_2023



# numDF denDF   F-value p-value
# (Intercept)     1    58 158.02196  <.0001
# TreatmentSB     3    58   2.46117  0.0716


#### Sample Graph For Legend ####
legend_2 <- ggplot(average_richness_2021, aes(x = mean_richness, color = "PBG")) +
  geom_density(alpha = 0.5) +
  geom_vline(aes(xintercept = mean_richness_ABG_2021 , color = "ABG"), linetype = "dashed", size = 1) +
  labs(title = "2021: Density Plot of Mean Richness",
       x = "Mean Richness",
       y = "Density") +
  scale_color_manual(values = c("blue", "red"), name = "Legend") +
  annotate("text", x = 13, y = 0.5, label = "p-value: 0.202", color = "blue", size = 3) +
  annotate("text", x = 13, y = 0.44, label = "Z-score: 1.276", color = "red", size = 3) +
  theme_bw()

legend <- get_legend(legend_2)

#### 2021-2023 Multi-panel Years Since Burned ####
# Arrange plots into a multipanel graph
theme_set(theme_bw())


multi_panel_graph <- grid.arrange(
  richness_boxplot, evenness_boxplot, count_boxplot,
  richness_boxplot_2022, evenness_boxplot_2022, count_boxplot_2022,
  richness_boxplot_2023, evenness_boxplot_2023, count_boxplot_2023,
  nrow = 3, ncol = 3)


# Display the multipanel graph
multi_panel_graph


#### Beta Diversity 2021 ####
set.seed(123)

library(vegan)
library(dplyr)


#Caculating Beta Diversity for PBG

Joined_New_2021 <- abundanceWide_2021 %>%
  filter(Treatment == "PBG") %>%
  group_by(block) %>%
  mutate(plot_index = row_number()) # Corrected index generation

num_bootstrap <- 1000
PBG_plot_master_2021 <- data.frame(iteration = 1:num_bootstrap, beta_diversity = numeric(num_bootstrap))

for (BOOT in 1:num_bootstrap) {
  Joined_New_Key <- Joined_New_2021 %>%
    sample_n(16, replace = TRUE) %>%
    select(plot_index) %>%
    ungroup()
  
  PBG_plot_ready <- Joined_New_2021 %>%
    right_join(Joined_New_Key, by = "plot_index") %>%
    mutate(iteration = BOOT)
  
  data_matrix <- PBG_plot_ready %>%
    select(8:152) # Species column
  
  dissimilarity_matrix <- vegdist(data_matrix, method = "bray")
  
  beta_diversity_value <- mean(dissimilarity_matrix)
  
  PBG_plot_master_2021$beta_diversity[BOOT] <- beta_diversity_value
}

#Add Treatment column

PBG_plot_master_2021$Treatment <- 'PBG'


# Print the results
print(PBG_plot_master_2021)

ggplot(PBG_plot_master_2021, aes(x = beta_diversity)) +
  geom_density() +
  labs(title = "Density Plot of Beta Diversity",
       x = "Beta Diversity",
       y = "Density") +
  theme_bw()



#Calculating ABG betadiveristy

Joined_New_2021_A <- abundanceWide_2021 %>%
  filter(Treatment == "ABG")

data_matrix_2 <- Joined_New_2021_A %>%
  select(8:151) # Species column

dissimilarity_matrix_2 <- vegdist(data_matrix, method = "bray")

beta_diversity_value_2 <- mean(dissimilarity_matrix)

#### Beta Diversity Graph 2021 ####


# Beta Diversity Density plot
ggplot(PBG_plot_master_2021, aes(x = beta_diversity)) +
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



#### Beta Diversity Z-score 2021 ####

#Getting average richness per iteration for bootstrapped dataframe
PBG_average_beta <- PBG_plot_master_2021 %>%
  summarize(mean_beta = mean(beta_diversity))


# Z-Score for richness 

Z_B <- ((beta_diversity_value_2) - (PBG_average_beta$mean_beta))/(sd(PBG_plot_master_2021$beta_diversity))
Z_B


p_value_B <- 1 - pnorm(Z_B)

#lower.tail = FALSE

print(p_value_B)

#Z-score = 0.2103348, p = 0.4167032





#### Beta Diversity 2022 ####
set.seed(123)

library(vegan)
library(dplyr)

# Calculating Beta Diversity for PBG

Joined_New_2022 <- abundanceWide_2022 %>%
  filter(Treatment == "PBG") %>%
  group_by(block) %>%
  mutate(plot_index = row_number()) # Corrected index generation

num_bootstrap <- 1000
PBG_plot_master_2022 <- data.frame(iteration = 1:num_bootstrap, beta_diversity = numeric(num_bootstrap))

for (BOOT in 1:num_bootstrap) {
  Joined_New_Key <- Joined_New_2022 %>%
    sample_n(16, replace = TRUE) %>%
    select(plot_index) %>%
    ungroup()
  
  PBG_plot_ready <- Joined_New_2022 %>%
    right_join(Joined_New_Key, by = "plot_index") %>%
    mutate(iteration = BOOT)
  
  data_matrix <- PBG_plot_ready %>%
    select(8:122) # Species column
  
  dissimilarity_matrix <- vegdist(data_matrix, method = "bray")
  
  beta_diversity_value <- mean(dissimilarity_matrix)
  
  PBG_plot_master_2022$beta_diversity[BOOT] <- beta_diversity_value
}

# Add Treatment column
PBG_plot_master_2022$Treatment <- 'PBG'

# Print the results
print(PBG_plot_master)

ggplot(PBG_plot_master_2022, aes(x = beta_diversity)) +
  geom_density() +
  labs(title = "Density Plot of Beta Diversity",
       x = "Beta Diversity",
       y = "Density") +
  theme_bw()

# Calculating ABG beta diversity

Joined_New_2022_A <- abundanceWide_2022 %>%
  filter(Treatment == "ABG")

data_matrix_2 <- Joined_New_2022_A %>%
  select(8:121) # Species column

dissimilarity_matrix_2 <- vegdist(data_matrix_2, method = "bray")

beta_diversity_value_2 <- mean(dissimilarity_matrix_2)

#### Beta Diversity Graph 2022 ####

# Beta Diversity Density plot
ggplot(PBG_plot_master_2022, aes(x = beta_diversity)) +
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

#### Beta Diversity Z-score 2022 ####

# Getting average richness per iteration for bootstrapped dataframe
PBG_average_beta <- PBG_plot_master_2022 %>%
  summarize(mean_beta = mean(beta_diversity))

# Z-Score for richness 

Z_B <- ((beta_diversity_value_2) - (PBG_average_beta$mean_beta))/(sd(PBG_plot_master_2022$beta_diversity))
Z_B

p_value_B <- 1 - pnorm(Z_B, lower.tail =  F)

# lower.tail = FALSE

print(p_value_B)

# Z-score = -0.01087649, p = 0.495661

#### Beta Diversity 2023 ####
set.seed(123)

library(vegan)
library(dplyr)

# Calculating Beta Diversity for PBG

Joined_New_2023 <- abundanceWide_2023 %>%
  filter(Treatment == "PBG") %>%
  group_by(block) %>%
  mutate(plot_index = row_number()) # Corrected index generation

num_bootstrap <- 1000
PBG_plot_master_2023 <- data.frame(iteration = 1:num_bootstrap, beta_diversity = numeric(num_bootstrap))

for (BOOT in 1:num_bootstrap) {
  Joined_New_Key <- Joined_New_2023 %>%
    sample_n(16, replace = TRUE) %>%
    select(plot_index) %>%
    ungroup()
  
  PBG_plot_ready <- Joined_New_2023 %>%
    right_join(Joined_New_Key, by = "plot_index") %>%
    mutate(iteration = BOOT)
  
  data_matrix <- PBG_plot_ready %>%
    select(8:77) # Species column
  
  dissimilarity_matrix <- vegdist(data_matrix, method = "bray")
  
  beta_diversity_value <- mean(dissimilarity_matrix)
  
  PBG_plot_master_2023$beta_diversity[BOOT] <- beta_diversity_value
}

# Add Treatment column
PBG_plot_master_2023$Treatment <- 'PBG'

# Print the results
print(PBG_plot_master_2023)

ggplot(PBG_plot_master_2023, aes(x = beta_diversity)) +
  geom_density() +
  labs(title = "Density Plot of Beta Diversity",
       x = "Beta Diversity",
       y = "Density") +
  theme_bw()

# Calculating ABG beta diversity

Joined_New_2023_A <- abundanceWide_2023 %>%
  filter(Treatment == "ABG")

data_matrix_2 <- Joined_New_2023_A %>%
  select(8:76) # Species column

dissimilarity_matrix_2 <- vegdist(data_matrix_2, method = "bray")

beta_diversity_value_2 <- mean(dissimilarity_matrix_2)

#### Beta Diversity Graph 2023 ####

# Beta Diversity Density plot
ggplot(PBG_plot_master_2023, aes(x = beta_diversity)) +
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

#### Beta Diversity Z-score 2023 ####

# Getting average richness per iteration for bootstrapped dataframe
PBG_average_beta <- PBG_plot_master_2023 %>%
  summarize(mean_beta = mean(beta_diversity))

# Z-Score for richness 

Z_B <- ((beta_diversity_value_2) - (PBG_average_beta$mean_beta))/(sd(PBG_plot_master_2023$beta_diversity))
Z_B

p_value_B <- 1 - pnorm(Z_B)

# lower.tail = FALSE

print(p_value_B)

# Z-score = 2.477193, p = 0.006621016
 

#### Functional Traits Unique Species Pull ####

#Grab unique species

# Get unique families
unique_families <- unique(Abundance_ID_aboveground$Family)

# Convert to data frame
unique_families_df <- data.frame(Family = unique_families)

# Write to Excel
write.xlsx(unique_families_df, file = "unique_families_aboveground.xlsx")

#### Functional Traits Cleanup Code ####
Functional_Groups <- read_excel("2021-2023 Aboveground PBG Inv Functional Groups.xlsx")

merged_data <- Abundance_ID_aboveground %>%
  left_join(Functional_Groups, by = c("Order", "Family"))

cleaned_data <- merged_data %>% filter(!is.na(Functional_Group))

#### Functional Trait Bootstrapping 2021 ####
# Filter data for 2021
data_2021 <- cleaned_data %>% filter(Date == 2021)

# Group by Sample and Functional_Group for 2021
abundance_summary_2021 <- data_2021 %>%
  group_by(Sample, Functional_Group) %>%
  summarise(Total_Abundance = sum(Count, na.rm = TRUE), .groups = 'drop')

# Add block and plot_index columns
abundance_summary_2021_block <- abundance_summary_2021 %>%
  mutate(
    block = ifelse(grepl("s", Sample, ignore.case = TRUE), "North", "South"),
    Treatment = case_when(
      grepl("1", Sample) ~ "ABG",
      grepl("3", Sample) ~ "PBG",
      TRUE ~ NA_character_
    )
  )

sample_wide_2021 <- abundance_summary_2021_block %>%
  pivot_wider(
    names_from = Functional_Group,
    values_from = Total_Abundance,
    values_fill = 0  # Fill missing values with 0
  )

sample_compact_2021 <- sample_wide_2021 %>%
  filter(Treatment == "PBG") %>% 
  group_by(Sample, block, Treatment) %>%
  summarise(across(Detritivore:`Gall Former`, sum, na.rm = TRUE), .groups = 'drop') %>%
  mutate(plot_index = 1:length(block))


num_bootstrap <- 1000
bootstrap_vector <- 1:num_bootstrap
Plot_master_2021<- data.frame()  # Initialize an empty dataframe


for (BOOT in bootstrap_vector) {
  # Sample unique plot_index values within each block for 2022
  New_Key_2021 <- sample_compact_2021 %>%
    dplyr::select(1:14) %>%
    unique() %>%
    group_by(block) %>%
    sample_n(16, replace = TRUE) %>%
    dplyr::select(plot_index) %>%
    ungroup()
  
  # Join the sampled rows back to the original dataframe
  Plot_ready_2021 <- sample_compact_2021 %>%
    right_join(New_Key_2021, by = c("block", "plot_index")) %>%
    mutate(iteration = BOOT)
  
  # Append the results to the master dataframe for 2022
  Plot_master_2021 <- rbind(Plot_master_2021, Plot_ready_2021)
}
#### Functional Trait Bootstrapping 2022 ####

# Filter data for 2022
data_2022 <- cleaned_data %>% filter(Date == 2022)

# Group by Sample and Functional_Group for 2022
abundance_summary_2022 <- data_2022 %>%
  group_by(Sample, Functional_Group) %>%
  summarise(Total_Abundance = sum(Count, na.rm = TRUE), .groups = 'drop')

# Add block and Treatment columns
abundance_summary_2022_block <- abundance_summary_2022 %>%
  mutate(
    block = ifelse(grepl("s", Sample, ignore.case = TRUE), "North", "South"),
    Treatment = case_when(
      grepl("1", Sample) ~ "ABG",
      grepl("3", Sample) ~ "PBG",
      TRUE ~ NA_character_
    )
  )

# Pivot to wide format based on Functional_Group and fill missing values with 0
sample_wide_2022 <- abundance_summary_2022_block %>%
  pivot_wider(
    names_from = Functional_Group,
    values_from = Total_Abundance,
    values_fill = 0
  )

# Filter for PBG treatment and calculate sum across functional groups
sample_compact_2022 <- sample_wide_2022 %>%
  filter(Treatment == "PBG") %>% 
  group_by(Sample, block, Treatment) %>%
  summarise(across(Fungivore:Parasite, sum, na.rm = TRUE), .groups = 'drop') %>%
  mutate(plot_index = 1:length(block))


num_bootstrap <- 1000
bootstrap_vector <- 1:num_bootstrap
Plot_master_2022 <- data.frame()  # Initialize an empty dataframe

for (BOOT in bootstrap_vector) {
  # Sample unique plot_index values within each block for 2022
  New_Key_2022 <- sample_compact_2022 %>%
    dplyr::select(1:12) %>%
    unique() %>%
    group_by(block) %>%
    sample_n(16, replace = TRUE) %>%
    dplyr::select(plot_index) %>%
    ungroup()
  
  # Join the sampled rows back to the original dataframe
  Plot_ready_2022 <- sample_compact_2022 %>%
    right_join(New_Key_2022, by = c("block", "plot_index")) %>%
    mutate(iteration = BOOT)
  
  # Append the results to the master dataframe for 2022
  Plot_master_2022 <- rbind(Plot_master_2022, Plot_ready_2022)
}

#### Functional Trait Bootstrapping 2023 ####

# Filter data for 2023
data_2023 <- cleaned_data %>% filter(Date == 2023)

# Group by Sample and Functional_Group for 2023
abundance_summary_2023 <- data_2023 %>%
  group_by(Sample, Functional_Group) %>%
  summarise(Total_Abundance = sum(Count, na.rm = TRUE), .groups = 'drop')

# Add block and Treatment columns
abundance_summary_2023_block <- abundance_summary_2022 %>%
  mutate(
    block = ifelse(grepl("s", Sample, ignore.case = TRUE), "North", "South"),
    Treatment = case_when(
      grepl("1", Sample) ~ "ABG",
      grepl("3", Sample) ~ "PBG",
      TRUE ~ NA_character_
    )
  )

# Pivot to wide format based on Functional_Group and fill missing values with 0
sample_wide_2023 <- abundance_summary_2022_block %>%
  pivot_wider(
    names_from = Functional_Group,
    values_from = Total_Abundance,
    values_fill = 0
  )

# Filter for PBG treatment and calculate sum across functional groups
sample_compact_2023 <- sample_wide_2022 %>%
  filter(Treatment == "PBG") %>% 
  group_by(Sample, block, Treatment) %>%
  summarise(across(Fungivore:Parasite, sum, na.rm = TRUE), .groups = 'drop') %>%
  mutate(plot_index = 1:length(block))

num_bootstrap <- 1000
bootstrap_vector <- 1:num_bootstrap
Plot_master_2023 <- data.frame()  # Initialize an empty dataframe

for (BOOT in bootstrap_vector) {
  # Sample unique plot_index values within each block for 2023
  New_Key_2023 <- sample_compact_2023 %>%
    dplyr::select(1:12) %>%
    unique() %>%
    group_by(block) %>%
    sample_n(16, replace = TRUE) %>%
    dplyr::select(plot_index) %>%
    ungroup()
  
  # Join the sampled rows back to the original dataframe
  Plot_ready_2023 <- sample_compact_2023 %>%
    right_join(New_Key_2023, by = c("block", "plot_index")) %>%
    mutate(iteration = BOOT)
  
  # Append the results to the master dataframe for 2023
  Plot_master_2023 <- rbind(Plot_master_2023, Plot_ready_2023)
}


#### 2021 Z-Score Calculations ####
detritivore_average_total_count_2021 <- Plot_master_2021 %>%
  group_by(iteration) %>%
  summarize(mean_count_detritivore = mean(Detritivore, na.rm = TRUE))

detritivore_mean_mean_total_count_2021 <- mean(detritivore_average_total_count_2021$mean_count_detritivore, na.rm = TRUE)

D_total_counts_2021 <- sample_wide_2021 %>%
  filter(str_ends(Sample, "2") | str_ends(Sample, "4")) %>%
  group_by(Sample, block, Treatment) %>%
  summarise(Count = sum(Detritivore, na.rm = TRUE), .groups = 'drop')

D_total_counts_ABG_2021 <- D_total_counts_2021 %>% filter(Treatment == "ABG")

D_ABG_mean_total_count_2021 <- mean(D_total_counts_ABG_2021$Count, na.rm = TRUE)


# Z-Score for total count 
D_C_2021 <- ((D_ABG_mean_total_count_2021) - (detritivore_mean_mean_total_count_2021))/(sd(detritivore_average_total_count_2021$mean_count_detritivore))
D_C_2021

p_value_D_C_2021 <- 2*pnorm(-abs(D_C_2021))
p_value_D_C_2021

#Z = -1.550103
#P = 0.1211169

# Herbivore
herbivore_average_total_count_2021 <- Plot_master_2021 %>%
  group_by(iteration) %>%
  summarize(mean_count_herbivore = mean(Herbivore, na.rm = TRUE))

herbivore_mean_mean_total_count_2021 <- mean(herbivore_average_total_count_2021$mean_count_herbivore, na.rm = TRUE)

D_total_counts_herbivore_2021 <- sample_wide_2021 %>%
  filter(str_ends(Sample, "2") | str_ends(Sample, "4")) %>%
  group_by(Sample, block, Treatment) %>%
  summarise(Count = sum(Herbivore, na.rm = TRUE), .groups = 'drop')

D_total_counts_herbivore_ABG_2021 <- D_total_counts_herbivore_2021 %>% filter(Treatment == "ABG")

D_ABG_mean_total_count_herbivore_2021 <- mean(D_total_counts_herbivore_ABG_2021$Count, na.rm = TRUE)

# Z-Score and p-value for Herbivore
D_C_herbivore_2021 <- (D_ABG_mean_total_count_herbivore_2021 - herbivore_mean_mean_total_count_2021) /
  sd(herbivore_average_total_count_2021$mean_count_herbivore)
D_C_herbivore_2021

p_value_D_C_herbivore_2021 <- 2 * pnorm(-abs(D_C_herbivore_2021))
p_value_D_C_herbivore_2021

#Z = 1.375491
#P =  0.1689794

# Omnivore
omnivore_average_total_count_2021 <- Plot_master_2021 %>%
  group_by(iteration) %>%
  summarize(mean_count_omnivore = mean(Omnivore, na.rm = TRUE))

omnivore_mean_mean_total_count_2021 <- mean(omnivore_average_total_count_2021$mean_count_omnivore, na.rm = TRUE)

D_total_counts_omnivore_2021 <- sample_wide_2021 %>%
  filter(str_ends(Sample, "2") | str_ends(Sample, "4")) %>%
  group_by(Sample, block, Treatment) %>%
  summarise(Count = sum(Omnivore, na.rm = TRUE), .groups = 'drop')

D_total_counts_omnivore_ABG_2021 <- D_total_counts_omnivore_2021 %>% filter(Treatment == "ABG")

D_ABG_mean_total_count_omnivore_2021 <- mean(D_total_counts_omnivore_ABG_2021$Count, na.rm = TRUE)

# Z-Score and p-value for Omnivore
D_C_omnivore_2021 <- (D_ABG_mean_total_count_omnivore_2021 - omnivore_mean_mean_total_count_2021) /
  sd(omnivore_average_total_count_2021$mean_count_omnivore)
D_C_omnivore_2021

p_value_D_C_omnivore_2021 <- 2 * pnorm(-abs(D_C_omnivore_2021))
p_value_D_C_omnivore_2021

#Z = 1.966422
#p = 0.04924985 *ABG Higher

# Predator
predator_average_total_count_2021 <- Plot_master_2021 %>%
  group_by(iteration) %>%
  summarize(mean_count_predator = mean(Predator, na.rm = TRUE))

predator_mean_mean_total_count_2021 <- mean(predator_average_total_count_2021$mean_count_predator, na.rm = TRUE)

D_total_counts_predator_2021 <- sample_wide_2021 %>%
  filter(str_ends(Sample, "2") | str_ends(Sample, "4")) %>%
  group_by(Sample, block, Treatment) %>%
  summarise(Count = sum(Predator, na.rm = TRUE), .groups = 'drop')

D_total_counts_predator_ABG_2021 <- D_total_counts_predator_2021 %>% filter(Treatment == "ABG")

D_ABG_mean_total_count_predator_2021 <- mean(D_total_counts_predator_ABG_2021$Count, na.rm = TRUE)

# Z-Score and p-value for Predator
D_C_predator_2021 <- (D_ABG_mean_total_count_predator_2021 - predator_mean_mean_total_count_2021) /
  sd(predator_average_total_count_2021$mean_count_predator)
D_C_predator_2021

p_value_D_C_predator_2021 <- 2 * pnorm(-abs(D_C_predator_2021))
p_value_D_C_predator_2021

# Z = 1.407006
# P = 0.1594256

# Fungivore
fungivore_average_total_count_2021 <- Plot_master_2021 %>%
  group_by(iteration) %>%
  summarize(mean_count_fungivore = mean(Fungivore, na.rm = TRUE))

fungivore_mean_mean_total_count_2021 <- mean(fungivore_average_total_count_2021$mean_count_fungivore, na.rm = TRUE)

D_total_counts_fungivore_2021 <- sample_wide_2021 %>%
  filter(str_ends(Sample, "2") | str_ends(Sample, "4")) %>%
  group_by(Sample, block, Treatment) %>%
  summarise(Count = sum(Fungivore, na.rm = TRUE), .groups = 'drop')

D_total_counts_fungivore_ABG_2021 <- D_total_counts_fungivore_2021 %>% filter(Treatment == "ABG")

D_ABG_mean_total_count_fungivore_2021 <- mean(D_total_counts_fungivore_ABG_2021$Count, na.rm = TRUE)

# Z-Score and p-value for Fungivore
D_C_fungivore_2021 <- (D_ABG_mean_total_count_fungivore_2021 - fungivore_mean_mean_total_count_2021) /
  sd(fungivore_average_total_count_2021$mean_count_fungivore)
D_C_fungivore_2021

p_value_D_C_fungivore_2021 <- 2 * pnorm(-abs(D_C_fungivore_2021))
p_value_D_C_fungivore_2021

# Z = -0.03054093
# P = 0.9756356

# Pollinator
pollinator_average_total_count_2021 <- Plot_master_2021 %>%
  group_by(iteration) %>%
  summarize(mean_count_pollinator = mean(Pollinator, na.rm = TRUE))

pollinator_mean_mean_total_count_2021 <- mean(pollinator_average_total_count_2021$mean_count_pollinator, na.rm = TRUE)

D_total_counts_pollinator_2021 <- sample_wide_2021 %>%
  filter(str_ends(Sample, "2") | str_ends(Sample, "4")) %>%
  group_by(Sample, block, Treatment) %>%
  summarise(Count = sum(Pollinator, na.rm = TRUE), .groups = 'drop')

D_total_counts_pollinator_ABG_2021 <- D_total_counts_pollinator_2021 %>% filter(Treatment == "ABG")

D_ABG_mean_total_count_pollinator_2021 <- mean(D_total_counts_pollinator_ABG_2021$Count, na.rm = TRUE)

# Z-Score and p-value for Pollinator
D_C_pollinator_2021 <- (D_ABG_mean_total_count_pollinator_2021 - pollinator_mean_mean_total_count_2021) /
  sd(pollinator_average_total_count_2021$mean_count_pollinator)
D_C_pollinator_2021

p_value_D_C_pollinator_2021 <- 2 * pnorm(-abs(D_C_pollinator_2021))
p_value_D_C_pollinator_2021

# Z = -1.566795
# P = 0.1171625

# Parasitoid
parasitoid_average_total_count_2021 <- Plot_master_2021 %>%
  group_by(iteration) %>%
  summarize(mean_count_parasitoid = mean(Parasitoid, na.rm = TRUE))

parasitoid_mean_mean_total_count_2021 <- mean(parasitoid_average_total_count_2021$mean_count_parasitoid, na.rm = TRUE)

D_total_counts_parasitoid_2021 <- sample_wide_2021 %>%
  group_by(Sample, block, Treatment) %>%
  summarise(Count = sum(Parasitoid, na.rm = TRUE), .groups = 'drop')

D_total_counts_parasitoid_ABG_2021 <- D_total_counts_parasitoid_2021 %>% filter(Treatment == "ABG")

D_ABG_mean_total_count_parasitoid_2021 <- mean(D_total_counts_parasitoid_ABG_2021$Count, na.rm = TRUE)

# Z-Score and p-value for Parasitoid
D_C_parasitoid_2021 <- (D_ABG_mean_total_count_parasitoid_2021 - parasitoid_mean_mean_total_count_2021) /
  sd(parasitoid_average_total_count_2021$mean_count_parasitoid)
D_C_parasitoid_2021

p_value_D_C_parasitoid_2021 <- 2 * pnorm(-abs(D_C_parasitoid_2021))
p_value_D_C_parasitoid_2021

# Z = -1.019089
# p = 0.3081605


# Parasite
parasite_average_total_count_2021 <- Plot_master_2021 %>%
  group_by(iteration) %>%
  summarize(mean_count_parasite = mean(Parasite, na.rm = TRUE))

parasite_mean_mean_total_count_2021 <- mean(parasite_average_total_count_2021$mean_count_parasite, na.rm = TRUE)

D_total_counts_parasite_2021 <- sample_wide_2021 %>%
  group_by(Sample, block, Treatment) %>%
  summarise(Count = sum(Parasite, na.rm = TRUE), .groups = 'drop')

D_total_counts_parasite_ABG_2021 <- D_total_counts_parasite_2021 %>% filter(Treatment == "ABG")

D_ABG_mean_total_count_parasite_2021 <- mean(D_total_counts_parasite_ABG_2021$Count, na.rm = TRUE)

# Z-Score and p-value for Parasite
D_C_parasite_2021 <- (D_ABG_mean_total_count_parasite_2021 - parasite_mean_mean_total_count_2021) /
  sd(parasite_average_total_count_2021$mean_count_parasite)
D_C_parasite_2021

p_value_D_C_parasite_2021 <- 2 * pnorm(-abs(D_C_parasite_2021))
p_value_D_C_parasite_2021

#Z = -1.439339
#P = 0.1500544


# Saprophagic
saprophagic_average_total_count_2021 <- Plot_master_2021 %>%
  group_by(iteration) %>%
  summarize(mean_count_saprophagic = mean(Saprophagic, na.rm = TRUE))

saprophagic_mean_mean_total_count_2021 <- mean(saprophagic_average_total_count_2021$mean_count_saprophagic, na.rm = TRUE)

D_total_counts_saprophagic_2021 <- sample_wide_2021 %>%
  group_by(Sample, block, Treatment) %>%
  summarise(Count = sum(Saprophagic, na.rm = TRUE), .groups = 'drop')

D_total_counts_saprophagic_ABG_2021 <- D_total_counts_saprophagic_2021 %>% filter(Treatment == "ABG")

D_ABG_mean_total_count_saprophagic_2021 <- mean(D_total_counts_saprophagic_ABG_2021$Count, na.rm = TRUE)

# Z-Score and p-value for Saprophagic
D_C_saprophagic_2021 <- (D_ABG_mean_total_count_saprophagic_2021 - saprophagic_mean_mean_total_count_2021) /
  sd(saprophagic_average_total_count_2021$mean_count_saprophagic)
D_C_saprophagic_2021

p_value_D_C_saprophagic_2021 <- 2 * pnorm(-abs(D_C_saprophagic_2021))
p_value_D_C_saprophagic_2021

# Z = -1.17403
# P = 0.2403829

# Gall Former
gall_former_average_total_count_2021 <- Plot_master_2021 %>%
  group_by(iteration) %>%
  summarize(mean_count_gall_former = mean(`Gall Former`, na.rm = TRUE))

gall_former_mean_mean_total_count_2021 <- mean(gall_former_average_total_count_2021$mean_count_gall_former, na.rm = TRUE)

D_total_counts_gall_former_2021 <- sample_wide_2021 %>%
  group_by(Sample, block, Treatment) %>%
  summarise(Count = sum(`Gall Former`, na.rm = TRUE), .groups = 'drop')

D_total_counts_gall_former_ABG_2021 <- D_total_counts_gall_former_2021 %>% filter(Treatment == "ABG")

D_ABG_mean_total_count_gall_former_2021 <- mean(D_total_counts_gall_former_ABG_2021$Count, na.rm = TRUE)

# Z-Score and p-value for Gall Former
D_C_gall_former_2021 <- (D_ABG_mean_total_count_gall_former_2021 - gall_former_mean_mean_total_count_2021) /
  sd(gall_former_average_total_count_2021$mean_count_gall_former)
D_C_gall_former_2021

p_value_D_C_gall_former_2021 <- 2 * pnorm(-abs(D_C_gall_former_2021))
p_value_D_C_gall_former_2021

#Z = -0.8498405
# P = 0.3954138

#### 2022 Z-score Calculations ####

# Detritivore
detritivore_average_total_count_2022 <- Plot_master_2022 %>%
  group_by(iteration) %>%
  summarize(mean_count_detritivore = mean(Detritivore, na.rm = TRUE))

detritivore_mean_mean_total_count_2022 <- mean(detritivore_average_total_count_2022$mean_count_detritivore, na.rm = TRUE)

D_total_counts_detritivore_2022 <- sample_wide_2022 %>%
  group_by(Sample, block, Treatment) %>%
  summarise(Count = sum(Detritivore, na.rm = TRUE), .groups = 'drop')

D_total_counts_detritivore_ABG_2022 <- D_total_counts_detritivore_2022 %>% filter(Treatment == "ABG")

D_ABG_mean_total_count_detritivore_2022 <- mean(D_total_counts_detritivore_ABG_2022$Count, na.rm = TRUE)

# Z-Score and p-value for Detritivore
D_C_detritivore_2022 <- (D_ABG_mean_total_count_detritivore_2022 - detritivore_mean_mean_total_count_2022) /
  sd(detritivore_average_total_count_2022$mean_count_detritivore)
D_C_detritivore_2022
p_value_D_C_detritivore_2022 <- 2 * pnorm(-abs(D_C_detritivore_2022))
p_value_D_C_detritivore_2022

#Z = -0.4957635
#P = 0.6200612

# Herbivore
herbivore_average_total_count_2022 <- Plot_master_2022 %>%
  group_by(iteration) %>%
  summarize(mean_count_herbivore = mean(Herbivore, na.rm = TRUE))

herbivore_mean_mean_total_count_2022 <- mean(herbivore_average_total_count_2022$mean_count_herbivore, na.rm = TRUE)

D_total_counts_herbivore_2022 <- sample_wide_2022 %>%
  group_by(Sample, block, Treatment) %>%
  summarise(Count = sum(Herbivore, na.rm = TRUE), .groups = 'drop')

D_total_counts_herbivore_ABG_2022 <- D_total_counts_herbivore_2022 %>% filter(Treatment == "ABG")

D_ABG_mean_total_count_herbivore_2022 <- mean(D_total_counts_herbivore_ABG_2022$Count, na.rm = TRUE)

# Z-Score and p-value for Herbivore
D_C_herbivore_2022 <- (D_ABG_mean_total_count_herbivore_2022 - herbivore_mean_mean_total_count_2022) /
  sd(herbivore_average_total_count_2022$mean_count_herbivore)
D_C_herbivore_2022

p_value_D_C_herbivore_2022 <- 2 * pnorm(-abs(D_C_herbivore_2022))
p_value_D_C_herbivore_2022

#Z = 1.011356
#P = 0.3118459

# Omnivore
omnivore_average_total_count_2022 <- Plot_master_2022 %>%
  group_by(iteration) %>%
  summarize(mean_count_omnivore = mean(Omnivore, na.rm = TRUE))

omnivore_mean_mean_total_count_2022 <- mean(omnivore_average_total_count_2022$mean_count_omnivore, na.rm = TRUE)

D_total_counts_omnivore_2022 <- sample_wide_2022 %>%
  group_by(Sample, block, Treatment) %>%
  summarise(Count = sum(Omnivore, na.rm = TRUE), .groups = 'drop')

D_total_counts_omnivore_ABG_2022 <- D_total_counts_omnivore_2022 %>% filter(Treatment == "ABG")

D_ABG_mean_total_count_omnivore_2022 <- mean(D_total_counts_omnivore_ABG_2022$Count, na.rm = TRUE)

# Z-Score and p-value for Omnivore
D_C_omnivore_2022 <- (D_ABG_mean_total_count_omnivore_2022 - omnivore_mean_mean_total_count_2022) /
  sd(omnivore_average_total_count_2022$mean_count_omnivore)
D_C_omnivore_2022

p_value_D_C_omnivore_2022 <- 2 * pnorm(-abs(D_C_omnivore_2022))
p_value_D_C_omnivore_2022

# Z = 0.7971701
# P = 0.4253522

# Predator
predator_average_total_count_2022 <- Plot_master_2022 %>%
  group_by(iteration) %>%
  summarize(mean_count_predator = mean(Predator, na.rm = TRUE))

predator_mean_mean_total_count_2022 <- mean(predator_average_total_count_2022$mean_count_predator, na.rm = TRUE)

D_total_counts_predator_2022 <- sample_wide_2022 %>%
  group_by(Sample, block, Treatment) %>%
  summarise(Count = sum(Predator, na.rm = TRUE), .groups = 'drop')

D_total_counts_predator_ABG_2022 <- D_total_counts_predator_2022 %>% filter(Treatment == "ABG")

D_ABG_mean_total_count_predator_2022 <- mean(D_total_counts_predator_ABG_2022$Count, na.rm = TRUE)

# Z-Score and p-value for Predator
D_C_predator_2022 <- (D_ABG_mean_total_count_predator_2022 - predator_mean_mean_total_count_2022) /
  sd(predator_average_total_count_2022$mean_count_predator)
D_C_predator_2022

p_value_D_C_predator_2022 <- 2 * pnorm(-abs(D_C_predator_2022))
p_value_D_C_predator_2022

#Z = 1.089409
#P = 0.2759736s

# Fungivore
fungivore_average_total_count_2022 <- Plot_master_2022 %>%
  group_by(iteration) %>%
  summarize(mean_count_fungivore = mean(Fungivore, na.rm = TRUE))

fungivore_mean_mean_total_count_2022 <- mean(fungivore_average_total_count_2022$mean_count_fungivore, na.rm = TRUE)

D_total_counts_fungivore_2022 <- sample_wide_2022 %>%
  group_by(Sample, block, Treatment) %>%
  summarise(Count = sum(Fungivore, na.rm = TRUE), .groups = 'drop')

D_total_counts_fungivore_ABG_2022 <- D_total_counts_fungivore_2022 %>% filter(Treatment == "ABG")

D_ABG_mean_total_count_fungivore_2022 <- mean(D_total_counts_fungivore_ABG_2022$Count, na.rm = TRUE)

# Z-Score and p-value for Fungivore
D_C_fungivore_2022 <- (D_ABG_mean_total_count_fungivore_2022 - fungivore_mean_mean_total_count_2022) /
  sd(fungivore_average_total_count_2022$mean_count_fungivore)
D_C_fungivore_2022

p_value_D_C_fungivore_2022 <- 2 * pnorm(-abs(D_C_fungivore_2022))
p_value_D_C_fungivore_2022

# Z = -1.273299
# P = 0.2029119

# Pollinator
pollinator_average_total_count_2022 <- Plot_master_2022 %>%
  group_by(iteration) %>%
  summarize(mean_count_pollinator = mean(Pollinator, na.rm = TRUE))

pollinator_mean_mean_total_count_2022 <- mean(pollinator_average_total_count_2022$mean_count_pollinator, na.rm = TRUE)

D_total_counts_pollinator_2022 <- sample_wide_2022 %>%
  group_by(Sample, block, Treatment) %>%
  summarise(Count = sum(Pollinator, na.rm = TRUE), .groups = 'drop')

D_total_counts_pollinator_ABG_2022 <- D_total_counts_pollinator_2022 %>% filter(Treatment == "ABG")

D_ABG_mean_total_count_pollinator_2022 <- mean(D_total_counts_pollinator_ABG_2022$Count, na.rm = TRUE)

# Z-Score and p-value for Pollinator
D_C_pollinator_2022 <- (D_ABG_mean_total_count_pollinator_2022 - pollinator_mean_mean_total_count_2022) /
  sd(pollinator_average_total_count_2022$mean_count_pollinator)
D_C_pollinator_2022

p_value_D_C_pollinator_2022 <- 2 * pnorm(-abs(D_C_pollinator_2022))
p_value_D_C_pollinator_2022

#Z = -1.307402
#P = 0.1910761

# Parasitoid
parasitoid_average_total_count_2022 <- Plot_master_2022 %>%
  group_by(iteration) %>%
  summarize(mean_count_parasitoid = mean(Parasitoid, na.rm = TRUE))

parasitoid_mean_mean_total_count_2022 <- mean(parasitoid_average_total_count_2022$mean_count_parasitoid, na.rm = TRUE)

D_total_counts_parasitoid_2022 <- sample_wide_2022 %>%
  group_by(Sample, block, Treatment) %>%
  summarise(Count = sum(Parasitoid, na.rm = TRUE), .groups = 'drop')

D_total_counts_parasitoid_ABG_2022 <- D_total_counts_parasitoid_2022 %>% filter(Treatment == "ABG")

D_ABG_mean_total_count_parasitoid_2022 <- mean(D_total_counts_parasitoid_ABG_2022$Count, na.rm = TRUE)

# Z-Score and p-value for Parasitoid
D_C_parasitoid_2022 <- (D_ABG_mean_total_count_parasitoid_2022 - parasitoid_mean_mean_total_count_2022) /
  sd(parasitoid_average_total_count_2022$mean_count_parasitoid)
D_C_parasitoid_2022

p_value_D_C_parasitoid_2022 <- 2 * pnorm(-abs(D_C_parasitoid_2022))
p_value_D_C_parasitoid_2022

# Z = -1.365106
# P = 0.1722199

# Parasite
parasite_average_total_count_2022 <- Plot_master_2022 %>%
  group_by(iteration) %>%
  summarize(mean_count_parasite = mean(Parasite, na.rm = TRUE))

parasite_mean_mean_total_count_2022 <- mean(parasite_average_total_count_2022$mean_count_parasite, na.rm = TRUE)

D_total_counts_parasite_2022 <- sample_wide_2022 %>%
  group_by(Sample, block, Treatment) %>%
  summarise(Count = sum(Parasite, na.rm = TRUE), .groups = 'drop')

D_total_counts_parasite_ABG_2022 <- D_total_counts_parasite_2022 %>% filter(Treatment == "ABG")

D_ABG_mean_total_count_parasite_2022 <- mean(D_total_counts_parasite_ABG_2022$Count, na.rm = TRUE)

# Z-Score and p-value for Parasite
D_C_parasite_2022 <- (D_ABG_mean_total_count_parasite_2022 - parasite_mean_mean_total_count_2022) /
  sd(parasite_average_total_count_2022$mean_count_parasite)
D_C_parasite_2022

p_value_D_C_parasite_2022 <- 2 * pnorm(-abs(D_C_parasite_2022))
p_value_D_C_parasite_2022

# Z = -0.8493654
# P = 0.395678
#### 2023 Z-score Calculations ####
# Detritivore
detritivore_average_total_count_2023 <- Plot_master_2023 %>%
  group_by(iteration) %>%
  summarize(mean_count_detritivore = mean(Detritivore, na.rm = TRUE))

detritivore_mean_mean_total_count_2023 <- mean(detritivore_average_total_count_2023$mean_count_detritivore, na.rm = TRUE)

D_total_counts_detritivore_2023 <- sample_wide_2023 %>%
  group_by(Sample, block,Treatment) %>%
  summarise(Count = sum(Detritivore, na.rm = TRUE), .groups = 'drop') %>% 
  filter(Treatment == "ABG")

D_ABG_mean_total_count_detritivore_2023 <- mean(D_total_counts_detritivore_2023$Count, na.rm = TRUE)


# Z-Score and p-value for Detritivore
D_C_detritivore_2023 <- (D_ABG_mean_total_count_detritivore_2023 - detritivore_mean_mean_total_count_2023) /
  sd(detritivore_average_total_count_2023$mean_count_detritivore)
D_C_detritivore_2023

p_value_D_C_detritivore_2023 <- 2 * pnorm(-abs(D_C_detritivore_2023))
p_value_D_C_detritivore_2023

# Z = -0.4758728
# P = 0.634165

# Herbivore
herbivore_average_total_count_2023 <- Plot_master_2023 %>%
  group_by(iteration) %>%
  summarize(mean_count_herbivore = mean(Herbivore, na.rm = TRUE))

herbivore_mean_mean_total_count_2023 <- mean(herbivore_average_total_count_2023$mean_count_herbivore, na.rm = TRUE)

D_total_counts_herbivore_2023 <- sample_wide_2023 %>%
  group_by(Sample, block, Treatment) %>%
  summarise(Count = sum(Herbivore, na.rm = TRUE), .groups = 'drop') %>% 
  filter(Treatment == "ABG")

D_ABG_mean_total_count_herbivore_2023 <- mean(D_total_counts_herbivore_2023$Count, na.rm = TRUE)

# Z-Score and p-value for Herbivore
D_C_herbivore_2023 <- (D_ABG_mean_total_count_herbivore_2023 - herbivore_mean_mean_total_count_2023) /
  sd(herbivore_average_total_count_2023$mean_count_herbivore)
D_C_herbivore_2023

p_value_D_C_herbivore_2023 <- 2 * pnorm(-abs(D_C_herbivore_2023))
p_value_D_C_herbivore_2023

#Z = 5.654136
#P = 1.566322e-08 *Higher in ABG

# Omnivore
omnivore_average_total_count_2023 <- Plot_master_2023 %>%
  group_by(iteration) %>%
  summarize(mean_count_omnivore = mean(Omnivore, na.rm = TRUE))

omnivore_mean_mean_total_count_2023 <- mean(omnivore_average_total_count_2023$mean_count_omnivore, na.rm = TRUE)

D_total_counts_omnivore_2023 <- sample_wide_2023 %>%
  group_by(Sample, block, Treatment) %>%
  summarise(Count = sum(Omnivore, na.rm = TRUE), .groups = 'drop') %>% 
  filter(Treatment == "ABG")

D_ABG_mean_total_count_omnivore_2023 <- mean(D_total_counts_omnivore_2023$Count, na.rm = TRUE)

# Z-Score and p-value for Omnivore
D_C_omnivore_2023 <- (D_ABG_mean_total_count_omnivore_2023 - omnivore_mean_mean_total_count_2023) /
  sd(omnivore_average_total_count_2023$mean_count_omnivore)
D_C_omnivore_2023

p_value_D_C_omnivore_2023 <- 2 * pnorm(-abs(D_C_omnivore_2023))
p_value_D_C_omnivore_2023

#Z = 0.7874125
#P = 0.4310404

# Predator
predator_average_total_count_2023 <- Plot_master_2023 %>%
  group_by(iteration) %>%
  summarize(mean_count_predator = mean(Predator, na.rm = TRUE))

predator_mean_mean_total_count_2023 <- mean(predator_average_total_count_2023$mean_count_predator, na.rm = TRUE)

D_total_counts_predator_2023 <- sample_wide_2023 %>%
  group_by(Sample, block, Treatment) %>%
  summarise(Count = sum(Predator, na.rm = TRUE), .groups = 'drop') %>% 
  filter(Treatment == "ABG")

D_ABG_mean_total_count_predator_2023 <- mean(D_total_counts_predator_2023$Count, na.rm = TRUE)

# Z-Score and p-value for Predator
D_C_predator_2023 <- (D_ABG_mean_total_count_predator_2023 - predator_mean_mean_total_count_2023) /
  sd(predator_average_total_count_2023$mean_count_predator)
D_C_predator_2023

p_value_D_C_predator_2023 <- 2 * pnorm(-abs(D_C_predator_2023))
p_value_D_C_predator_2023

# Z = 1.144022
# P = 0.2526144

# Fungivore
fungivore_average_total_count_2023 <- Plot_master_2023 %>%
  group_by(iteration) %>%
  summarize(mean_count_fungivore = mean(Fungivore, na.rm = TRUE))

fungivore_mean_mean_total_count_2023 <- mean(fungivore_average_total_count_2023$mean_count_fungivore, na.rm = TRUE)

D_total_counts_fungivore_2023 <- sample_wide_2023 %>%
  group_by(Sample, block, Treatment) %>%
  summarise(Count = sum(Fungivore, na.rm = TRUE), .groups = 'drop') %>% 
  filter(Treatment == "ABG")

D_ABG_mean_total_count_fungivore_2023 <- mean(D_total_counts_fungivore_2023$Count, na.rm = TRUE)

# Z-Score and p-value for Fungivore
D_C_fungivore_2023 <- (D_ABG_mean_total_count_fungivore_2023 - fungivore_mean_mean_total_count_2023) /
  sd(fungivore_average_total_count_2023$mean_count_fungivore)
D_C_fungivore_2023

p_value_D_C_fungivore_2023 <- 2 * pnorm(-abs(D_C_fungivore_2023))
p_value_D_C_fungivore_2023

# Z = -1.277524
# P = 0.2014173

# Pollinator
pollinator_average_total_count_2023 <- Plot_master_2023 %>%
  group_by(iteration) %>%
  summarize(mean_count_pollinator = mean(Pollinator, na.rm = TRUE))

pollinator_mean_mean_total_count_2023 <- mean(pollinator_average_total_count_2023$mean_count_pollinator, na.rm = TRUE)

D_total_counts_pollinator_2023 <- sample_wide_2023 %>%
  group_by(Sample, block, Treatment) %>%
  summarise(Count = sum(Pollinator, na.rm = TRUE), .groups = 'drop') %>% 
  filter(Treatment == "ABG")

D_ABG_mean_total_count_pollinator_2023 <- mean(D_total_counts_pollinator_2023$Count, na.rm = TRUE)

# Z-Score and p-value for Pollinator
D_C_pollinator_2023 <- (D_ABG_mean_total_count_pollinator_2023 - pollinator_mean_mean_total_count_2023) /
  sd(pollinator_average_total_count_2023$mean_count_pollinator)
D_C_pollinator_2023

p_value_D_C_pollinator_2023 <- 2 * pnorm(-abs(D_C_pollinator_2023))
p_value_D_C_pollinator_2023

#Z = -1.275276
# P = 0.2022117

# Parasitoid
parasitoid_average_total_count_2023 <- Plot_master_2023 %>%
  group_by(iteration) %>%
  summarize(mean_count_parasitoid = mean(Parasitoid, na.rm = TRUE))

parasitoid_mean_mean_total_count_2023 <- mean(parasitoid_average_total_count_2023$mean_count_parasitoid, na.rm = TRUE)

D_total_counts_parasitoid_2023 <- sample_wide_2023 %>%
  group_by(Sample, block, Treatment) %>%
  summarise(Count = sum(Parasitoid, na.rm = TRUE), .groups = 'drop') %>% 
  filter(Treatment == "ABG")

D_ABG_mean_total_count_parasitoid_2023 <- mean(D_total_counts_parasitoid_2023$Count, na.rm = TRUE)

# Z-Score and p-value for Parasitoid
D_C_parasitoid_2023 <- (D_ABG_mean_total_count_parasitoid_2023 - parasitoid_mean_mean_total_count_2023) /
  sd(parasitoid_average_total_count_2023$mean_count_parasitoid)
D_C_parasitoid_2023

p_value_D_C_parasitoid_2023 <- 2 * pnorm(-abs(D_C_parasitoid_2023))
p_value_D_C_parasitoid_2023

# Z = -1.255643
# P = 0.2092454

# Parasite
parasite_average_total_count_2023 <- Plot_master_2023 %>%
  group_by(iteration) %>%
  summarize(mean_count_parasite = mean(Parasite, na.rm = TRUE))

parasite_mean_mean_total_count_2023 <- mean(parasite_average_total_count_2023$mean_count_parasite, na.rm = TRUE)

D_total_counts_parasite_2023 <- sample_wide_2023 %>%
  group_by(Sample, block, Treatment) %>%
  summarise(Count = sum(Parasite, na.rm = TRUE), .groups = 'drop') %>% 
  filter(Treatment == "ABG")

D_ABG_mean_total_count_parasite_2023 <- mean(D_total_counts_parasite_2023$Count, na.rm = TRUE)

# Z-Score and p-value for Parasite
D_C_parasite_2023 <- (D_ABG_mean_total_count_parasite_2023 - parasite_mean_mean_total_count_2023) /
  sd(parasite_average_total_count_2023$mean_count_parasite)
D_C_parasite_2023

p_value_D_C_parasite_2023 <- 2 * pnorm(-abs(D_C_parasite_2023))
p_value_D_C_parasite_2023

# Z = -0.8419705#
# P = 0.3998045

