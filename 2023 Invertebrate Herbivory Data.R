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

#### Data Cleanup ####

Invertebrate_Herbivory_2023 <- Invertebrate_Herbivory_2023 %>% mutate(Treatment = ifelse(grepl("C1", WS), "ABG", "PBG"),
       Block = ifelse(grepl("S", WS), "North", "South")) %>% 
  unite("Sample", c("WS", "Trans", "Rep"), sep = "_", remove = FALSE) %>% 
  unite("Plant", c("Gensus", "Species"), sep = "_", remove = FALSE) %>% 
  rename(Herbivory = 'Herbivory (%)' )
  



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