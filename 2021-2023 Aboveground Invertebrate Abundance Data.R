#### Libraries ####
library(readxl)
library(tidyverse)
library(codyn)
library(ggplot2)
library(ggpubr)
library(gridExtra)

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


#### Aboveground Weight Data Clean up ####

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




#### 2021 PERMANOVA & NMDS ####
### by Years since burn

#Combing count with burn info


#merge with burn data
Abundance_Data2021 <- full_join(Abundance_ID_aboveground, BurnInfo2021, by = "WS") %>% 
  unite("TreatmentSB",c("Treatment","SB"), sep="_") %>% 
  filter(Date == "2021")

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
  mutate(sum = rowSums(abundanceWide[, c(8:166)], na.rm = TRUE)) %>% 
  filter(sum > 0, !(Sample %in% c('C1A_A_38_ABG', 'C3SA_D_16_PBG', 'C3SA_C_38_PBG')))


print(permanova <- adonis2(formula = abundanceWide[,8:166]~TreatmentSB, data=abundanceWide, permutations=999, method="bray"))
#F=1.6053, df=3,52, p=0.016

#betadisper
veg <- vegdist(abundanceWide[,8:167], method = "bray")
dispersion <- betadisper(veg, abundanceWide$TreatmentSB)
permutest(dispersion, pairwise=TRUE, permutations=999) 
#F=1.5519, df=3,52, p=0.213

BC_Data <- metaMDS(abundanceWide[,8:167])
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

### by watershed
# PERMANOVA
print(permanova <- adonis2(formula = abundanceWide[,8:167]~Treatment, data=abundanceWide, permutations=999, method="bray"))
#F=1.4498, df=1,52, p=0.119

#betadisper
veg <- vegdist(abundanceWide[,8:167], method = "bray")
dispersion <- betadisper(veg, abundanceWide$Treatment)
permutest(dispersion, pairwise=TRUE, permutations=999) 
#F=9e-4, df=1,52, p=0.976

BC_Data <- metaMDS(abundanceWide[,8:167])
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
ggplot(BC_NMDS_Graph, aes(x = MDS1, y = MDS2, color = group, linetype = group, shape = group)) +
  geom_point(size = 6) + 
  geom_path(data = BC_Ellipses, aes(x = NMDS1, y = NMDS2), size = 3) +
  labs(color = "", linetype = "", shape = "") +
  scale_colour_manual(values = c("blue", "red"), name = "") +
  scale_linetype_manual(values = c("solid", "twodash"), name = "") +
  scale_shape_manual(values = c(19, 19)) +
  xlab("NMDS1") + 
  ylab("NMDS2") + 
  annotate('text', x = min(BC_NMDS_Graph$MDS1) - 0.04, y = max(BC_NMDS_Graph$MDS2) + 0.03, label = expression(paste('F'['3,52'], ' = 1.61')), size = 8, hjust = 'left') +
  annotate('text', x = min(BC_NMDS_Graph$MDS1) - 0.04, y = max(BC_NMDS_Graph$MDS2) - 0.04, label = 'p = 0.016', size = 8, hjust = 'left') +
  theme_classic() +
  theme(axis.text.x = element_text(size = 24, color = "black"), 
        axis.text.y = element_text(size = 24, color = "black"), 
        axis.title.x = element_text(size = 24, color = 'black'),
        axis.title.y = element_text(size = 24, color = 'black'),
        legend.text = element_text(size = 24),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()
  )


