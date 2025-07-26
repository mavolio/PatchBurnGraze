#### Libraries ####
library(readxl)
library(tidyverse)
library(lme4)
library(emmeans)
library(RVAideMemoire) #backtransform emmeans
library(glmmTMB) #mixed models for zero-inflated count data (bug data)
library(ggpubr)
library(RColorBrewer)
library(tidyverse)
library(car)
library(gridExtra)

# Add end a model for TreatmentSB (ABG + PBG 1-3) 


#### CSV read ####
Invertebrate_Herbivory_2023 <- read_excel("HerbivoryDataSheetPBG.xlsx")
BurnInfo2023 <- read_excel("YearsSenseBurned2023.xlsx")
#### BurnInfo Add In ####

Invertebrate_Herbivory_2023 <- left_join(Invertebrate_Herbivory_2023, BurnInfo2023, by = "WS")

#### Data Cleanup ####

Invertebrate_Herbivory_2023 <- Invertebrate_Herbivory_2023 %>% mutate(Treatment = ifelse(grepl("C1", WS), "ABG", "PBG"),
                                                                      Block = ifelse(grepl("S", WS), "North", "South")) %>% 
  unite("Sample", c("WS", "Trans", "Rep"), sep = "_", remove = FALSE) %>% 
  unite("Plant", c("Gensus", "Species"), sep = "_", remove = FALSE) %>% 
  rename(Herbivory = 'Herbivory (%)' ) %>% 
  unite("TreatmentSB",c("Treatment","SB"), sep="_") %>% 
  mutate(Treatment = ifelse(grepl("C1", WS), "ABG", "PBG")) 

Invertebrate_Herbivory_2023 <- Invertebrate_Herbivory_2023 %>%  
  mutate(Treatment = ifelse(grepl("C1", Sample), "ABG", "PBG"),
         block = ifelse(grepl("S", Sample), "North", "South"))

#### Bootstrapping ####
# num_bootstrap <- 1000
# 
# #Remove plant, species
# 
# Invertebrate_Herbivory_2023_New <- Invertebrate_Herbivory_2023 %>% 
#   select(-Plant, -Gensus, -Species, -Herbivory, -Notes) %>% 
#   unique() %>%
#   filter(Treatment == "PBG") %>% 
#   mutate(plot_index=1:length(block)) 
# 
# bootstrap_vector <- 1:num_bootstrap
# PBG_rep_master_2023 <- data.frame()  # Initialize an empty dataframe
# 
# for (BOOT in bootstrap_vector) {
#   Joined_New_2023 <- Invertebrate_Herbivory_2023_New %>%
#     dplyr::select(1:9) %>%
#     unique() %>%
#     group_by(Block) %>%
#     sample_n(24, replace = TRUE) %>%
#     dplyr::select(plot_index) %>%
#     ungroup()
#   
#   # Join the sampled rows back to the original dataframe
#   PBG_plot_ready_2023 <- Invertebrate_Herbivory_2023_New %>%
#     right_join(Joined_New_2023, by = c("Block", "plot_index")) %>%
#     mutate(iteration = BOOT)
#   
#   # Append the results to the master dataframe
#   PBG_rep_master_2023 <- rbind(PBG_rep_master_2023, PBG_plot_ready_2023)
# }


#Merge Output with original full join?




#### Condensing down ####
# PBG_rep_master_2023_new <- PBG_rep_master_2023 %>% 
#   group_by(iteration, Plant) %>% 
#   mutate(Seq = row_number()) %>% 
#   ungroup()
# 
# PBG_rep_master_2023_new_2 <- PBG_rep_master_2023_new %>%
#   group_by(Plant, Seq) %>%
#   summarise(Herbivory = mean(Herbivory))


############################################################
##### continuous response variable (ABG vs PBG) #####
############################################################

#Make it binary

# Define the threshold value
threshold <- 0.1  # Adjust this threshold as needed

# Convert 'Herbivory (%)' into a binary variable
Invertebrate_Herbivory_2023$Damage <- ifelse(Invertebrate_Herbivory_2023$Herbivory > threshold, 1, 0)

# Part 1: probability of zero damage

summary(zero_model <- glmer(`Damage` ~ Treatment + (1|Plant),
                            data = Invertebrate_Herbivory_2023,
                            family = binomial))


Anova(zero_model, type='III') 
back.emmeans(emmeans(zero_model, ~Treatment), transform='log')
# Part 2: distribution of the continuous, non-zero data

Damage_2023 <- Invertebrate_Herbivory_2023 %>% filter(Damage == 1)

summary(cont_model <- lmer(log(Herbivory) ~ Treatment + (1|Plant),
                           data = Damage_2023))
Anova(cont_model) 
back.emmeans(emmeans(cont_model, ~Treatment), transform='log')
#back.emmeans(emmeans(cont_model, ~Plant), transform='log')


## Figure ##
EBpalette1 <- c("#228833", "#725DEF", "#EE6677", "#702082")
EBpalettewoforest <- c("#DC267F", "#785EF0")

ggplot(data=Invertebrate_Herbivory_2023, aes(x=Treatment, y=Herbivory, fill=Plant, color=Plant)) + 
  geom_violin(aes(fill=Plant, color=Plant,
                  fill=after_scale(colorspace::lighten(fill, .3))),
              size=1, bw=.3)
############################################################
##### continuous response variable (ABG vs PBG years since burned) #####
############################################################

#Make it binary

# Define the threshold value
threshold <- 0.1  # Adjust this threshold as needed

# Convert 'Herbivory (%)' into a binary variable
Invertebrate_Herbivory_2023$Damage <- ifelse(Invertebrate_Herbivory_2023$Herbivory > threshold, 1, 0)

# Part 1: probability of zero damage

summary(zero_model <- glmer(`Damage` ~ TreatmentSB + (1|Plant),
                            data = Invertebrate_Herbivory_2023,
                            family = binomial))
Anova(zero_model, type='III') 
back.emmeans(emmeans(zero_model, ~TreatmentSB), transform='log')
pairwise_results <- emmeans(zero_model, pairwise ~ TreatmentSB, adjust = "tukey")
pairwise_results

# Part 2: distribution of the continuous, non-zero data

Damage_2023 <- Invertebrate_Herbivory_2023 %>% filter(Damage == 1)

summary(cont_model <- lmer(log(Herbivory) ~ TreatmentSB + (1|Plant),
                           data = Damage_2023))
Anova(cont_model) 
back.emmeans(emmeans(cont_model, ~TreatmentSB), transform='log')
#back.emmeans(emmeans(cont_model, ~Plant), transform='log')
pairwise_results <- emmeans(cont_model, pairwise ~ TreatmentSB, adjust = "tukey")
pairwise_results

## Figure ##

ggplot(data=Invertebrate_Herbivory_2023, aes(x=TreatmentSB, y=Herbivory, fill=Plant, color=Plant)) + 
  geom_violin(aes(fill=Plant, color=Plant,
                  fill=after_scale(colorspace::lighten(fill, .3))),
              size=1, bw=.3)


#### New Graphs ####

ABGvsPBG <- ggplot(data=Invertebrate_Herbivory_2023, aes(x=Treatment, y=Herbivory)) +  
  geom_boxplot(aes(fill=Treatment), width=0.1, position=position_dodge(width=0.9), outlier.shape=NA) + 
  geom_violin(aes(fill=Treatment), size=1, bw=.3, position=position_dodge(width=0.9)) +
  geom_jitter(position=position_jitter(width=0.2, height=0, seed=123), alpha=0.5) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(), 
    legend.position = "none",
    legend.text = element_text(size = 15),     # Increases legend text size
    legend.title = element_text(size = 0),    # Increases legend title size
    axis.title = element_text(size = 20),      # Increases axis titles size
    axis.text = element_text(size = 15)        # Increases axis text size
  ) +
  scale_fill_manual(values=c("ABG"="blue", "PBG"="red")) +
  labs(title="", x="", y="")
  

ABGvsPBG



YearSinceBurned <- ggplot(data=Invertebrate_Herbivory_2023, aes(x=TreatmentSB, y=Herbivory, fill=TreatmentSB)) + 
  geom_violin(color="black", size=1, bw=.3) +
  geom_boxplot(width=0.1, position=position_dodge(width=0.9), outlier.shape=NA) + 
  geom_jitter(position=position_jitter(width=0.2, height=0, seed=123), alpha=0.5) +
  theme_minimal() +
  scale_fill_manual(
    values=c("ABG_0"="blue", "PBG_0"="red", "PBG_1"="red", "PBG_2"="red"),
    labels=c("ABG_0"="ABG 0", "PBG_0"="PBG 0", "PBG_1"="PBG 1", "PBG_2"="PBG 2")  # Custom labels for the legend
  ) +
  scale_x_discrete(labels=c("ABG_0"="ABG 0", "PBG_0"="PBG 0", "PBG_1"="PBG 1", "PBG_2"="PBG 2")) +  # Custom labels for the x-axis
  labs(fill="Treatment") +   # Changes legend title to "Treatment"
  theme(
    panel.grid = element_blank(), 
    legend.position = "none",        # Moves legend to upper left corner
    legend.text = element_text(size = 15),     # Increases legend text size
    legend.title = element_text(size = 18),    # Increases legend title size
    axis.title = element_text(size = 20),      # Increases axis titles size
    axis.text = element_text(size = 15)        # Increases axis text (tick labels) size
  ) + 
  labs(x = "", y = "")

YearSinceBurned

#### Percent Damage ####
# Total number of rows
total_rows <- nrow(Invertebrate_Herbivory_2023)

# Number of rows with damage
damage_rows <- nrow(Invertebrate_Herbivory_2023 %>% filter(Damage > 0))

# Percentage of rows with damage
percentage_damage <- (damage_rows / total_rows) * 100

percentage_damage


# Calculate total damage and percent damage by species
damage_by_species <- Invertebrate_Herbivory_2023 %>%
  group_by(Species) %>%
  summarise(
    Total_Damage = sum(Damage, na.rm = TRUE),
    Herbivory_Records = n()
  ) %>%
  mutate(Percent_Damage = (Total_Damage / Herbivory_Records) * 100)

# View the results
print(damage_by_species)

#### Combined Graph ####

multi_panel_graph <- grid.arrange(ABGvsPBG, YearSinceBurned, nrow = 1)
#### New Graphs ####

# ABG vs PBG Plot
ABGvsPBG <- ggplot(data=Invertebrate_Herbivory_2023, aes(x=Treatment, y=Herbivory)) +  
  geom_boxplot(aes(fill=Treatment), width=0.1, position=position_dodge(width=0.9), outlier.shape=NA) + 
  geom_violin(aes(fill=Treatment), size=1, bw=.3, position=position_dodge(width=0.9)) +
  geom_jitter(position=position_jitter(width=0.2, height=0, seed=123), alpha=0.5) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(), 
    legend.position = "none",
    legend.text = element_text(size = 18),      
    legend.title = element_text(size = 0),      
   # axis.title = element_text(size = 24),       # Increased axis label size
    axis.text = element_text(size = 30)         # Increased tick label size
  ) +
  scale_fill_manual(values=c("ABG"="blue", "PBG"="red")) +
  labs(x="", 
  y="")  
# Ensuring clear axis labels


# Year Since Burned Plot
YearSinceBurned <- ggplot(data=Invertebrate_Herbivory_2023, aes(x=TreatmentSB, y=Herbivory, fill=TreatmentSB)) + 
  geom_violin(color="black", size=1, bw=.3) +
  geom_boxplot(width=0.1, position=position_dodge(width=0.9), outlier.shape=NA) + 
  geom_jitter(position=position_jitter(width=0.2, height=0, seed=123), alpha=0.5) +
  theme_minimal() +
  scale_fill_manual(
    values=c("ABG_0"="blue", "PBG_0"="red", "PBG_1"="red", "PBG_2"="red"),
    labels=c("ABG_0"="ABG 0", "PBG_0"="PBG 0", "PBG_1"="PBG 1", "PBG_2"="PBG 2")
  ) +
  scale_x_discrete(labels=c("ABG_0"="ABG 0", "PBG_0"="PBG 0", "PBG_1"="PBG 1", "PBG_2"="PBG 2")) +  
  labs(fill="Treatment") +  
  theme(
    panel.grid = element_blank(), 
    legend.position = "none",       
    legend.text = element_text(size = 18),      
    legend.title = element_text(size = 20),     
   # axis.title = element_text(size = 24),       # Increased axis label size
    axis.text = element_text(size = 30)         # Increased tick label size
  ) + 
  labs(x = "", 
     y = ""
       )

# Display graphs side by side
grid.arrange(ABGvsPBG, YearSinceBurned, nrow = 1)

#### New Graph for Manuscript ####
library(ggplot2)
library(dplyr)
library(gridExtra)

# Recalculate binary Damage
Invertebrate_Herbivory_2023 <- Invertebrate_Herbivory_2023 %>%
  mutate(Damage = ifelse(Herbivory > 0.1, 1, 0),
         TreatmentSB_clean = gsub("_", " ", TreatmentSB))  # Replace underscores

## ---- PANEL A (only legend, inset top-left, no box) ----
panel_A <- Invertebrate_Herbivory_2023 %>%
  group_by(Treatment, Sample) %>%
  summarise(Prop_Damage = mean(Damage), .groups = "drop") %>%
  ggplot(aes(x = Treatment, y = Prop_Damage, fill = Treatment)) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4) +
  scale_fill_manual(values = c("ABG" = "blue", "PBG" = "red")) +
  annotate("text", x = Inf, y = Inf, label = "A", hjust = 1.3, vjust = 1.5, size = 6) +
  labs(x = "", y = "Proportion with Damage", fill = "Treatment") +
  theme_minimal() +
  theme(
    legend.position = c(0.05, 0.95),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank(),
    axis.text = element_text(size = 12),
    title = element_text(size = 14)
  )

## ---- PANEL B (no legend, cleaned x-axis) ----
panel_B <- Invertebrate_Herbivory_2023 %>%
  group_by(TreatmentSB_clean, Sample) %>%
  summarise(Prop_Damage = mean(Damage), .groups = "drop") %>%
  ggplot(aes(x = TreatmentSB_clean, y = Prop_Damage, fill = TreatmentSB_clean)) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4) +
  scale_fill_manual(values = c("ABG 0"="blue", "PBG 0"="red", "PBG 1"="red", "PBG 2"="red")) +
  annotate("text", x = Inf, y = Inf, label = "B", hjust = 1.3, vjust = 1.5, size = 6) +
  labs(x = "", y = "Proportion with Damage") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        title = element_text(size = 14))

## ---- PANEL C (no legend) ----
panel_C <- ggplot(filter(Invertebrate_Herbivory_2023, Damage == 1),
                  aes(x = Treatment, y = Herbivory, fill = Treatment)) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4) +
  scale_fill_manual(values = c("ABG" = "blue", "PBG" = "red")) +
  annotate("text", x = Inf, y = Inf, label = "C", hjust = 1.3, vjust = 1.5, size = 6) +
  labs(x = "", y = "Leaf Damage (%)") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        title = element_text(size = 14))

## ---- PANEL D (no legend, cleaned x-axis) ----
panel_D <- ggplot(filter(Invertebrate_Herbivory_2023, Damage == 1),
                  aes(x = TreatmentSB_clean, y = Herbivory, fill = TreatmentSB_clean)) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4) +
  scale_fill_manual(values = c("ABG 0"="blue", "PBG 0"="red", "PBG 1"="red", "PBG 2"="red")) +
  annotate("text", x = Inf, y = Inf, label = "D", hjust = 1.3, vjust = 1.5, size = 6) +
  labs(x = "", y = "Leaf Damage (%)") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        title = element_text(size = 14))

## ---- Combine & Save ----
combined_figure <- grid.arrange(panel_A, panel_B, panel_C, panel_D, nrow = 2)
ggsave("Figure_Herbivory_Panels.png", plot = combined_figure, width = 10, height = 8, dpi = 300)
