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
num_bootstrap <- 1000

#Remove plant, species

Invertebrate_Herbivory_2023_New <- Invertebrate_Herbivory_2023 %>% 
  select(-Plant, -Gensus, -Species, -Herbivory, -Notes) %>% 
  unique() %>%
  filter(Treatment == "PBG") %>% 
  mutate(plot_index=1:length(block)) 

bootstrap_vector <- 1:num_bootstrap
PBG_rep_master_2023 <- data.frame()  # Initialize an empty dataframe

for (BOOT in bootstrap_vector) {
  Joined_New_2023 <- Invertebrate_Herbivory_2023_New %>%
    dplyr::select(1:9) %>%
    unique() %>%
    group_by(Block) %>%
    sample_n(24, replace = TRUE) %>%
    dplyr::select(plot_index) %>%
    ungroup()
  
  # Join the sampled rows back to the original dataframe
  PBG_plot_ready_2023 <- Invertebrate_Herbivory_2023_New %>%
    right_join(Joined_New_2023, by = c("Block", "plot_index")) %>%
    mutate(iteration = BOOT)
  
  # Append the results to the master dataframe
  PBG_rep_master_2023 <- rbind(PBG_rep_master_2023, PBG_plot_ready_2023)
}


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

summary(zero_model <- glmer(`Damage` ~ Treatment*Plant + (1|Block),
                            data = Invertebrate_Herbivory_2023,
                            family = binomial))


Anova(zero_model, type='III') 

# Part 2: distribution of the continuous, non-zero data

Damage_2023 <- Invertebrate_Herbivory_2023 %>% filter(Damage == 1)

summary(cont_model <- lmer(log(Herbivory) ~ Treatment*Plant + (1|Block),
                           data = Damage_2023))
Anova(cont_model) 
back.emmeans(emmeans(cont_model, ~Treatment), transform='log')
back.emmeans(emmeans(cont_model, ~Plant), transform='log')


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

summary(zero_model <- glmer(`Damage` ~ TreatmentSB*Plant + (1|Block),
                            data = Invertebrate_Herbivory_2023,
                            family = binomial))


Anova(zero_model, type='III') 

# Part 2: distribution of the continuous, non-zero data

Damage_2023 <- Invertebrate_Herbivory_2023 %>% filter(Damage == 1)

summary(cont_model <- lmer(log(Herbivory) ~ TreatmentSB*Plant + (1|Block),
                           data = Damage_2023))
Anova(cont_model) 
back.emmeans(emmeans(cont_model, ~Treatment), transform='log')
back.emmeans(emmeans(cont_model, ~Plant), transform='log')


## Figure ##

ggplot(data=Invertebrate_Herbivory_2023, aes(x=TreatmentSB, y=Herbivory, fill=Plant, color=Plant)) + 
  geom_violin(aes(fill=Plant, color=Plant,
                  fill=after_scale(colorspace::lighten(fill, .3))),
              size=1, bw=.3)
