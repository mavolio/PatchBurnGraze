####Libraries ####


library(ggplot2)

library(readxl)

library(dplyr)

####Cleanup ####

PBG_2021_Weight_Below <- read_excel("2021 Belowground Weight Data.xlsx")

table(PBG_2021_Weight_Below$WS)

#Clean up####


#First is blank sense we want all rows. This eliminates the empty rows with sheet comments. 


PBG_2021_Weight_Below <- PBG_2021_Weight_Below[, c("WS","Trans", "Dist.", "Weight (mg)")]

#Adds a new column to tag "1" with ABG and "3" with PBG. 

PBG_2021_Weight_Below$Treatment <- ifelse(grepl("1", PBG_2021_Weight_Below$WS), "ABG", "PBG")

#Clean up WS

PBG_2021_Weight_Below$WS <- gsub("C35A","C3SA", PBG_2021_Weight_Below$WS) 
PBG_2021_Weight_Below$WS <- gsub("C35B","C3SB", PBG_2021_Weight_Below$WS)
PBG_2021_Weight_Below$WS <- gsub("C3CB", "C3C", PBG_2021_Weight_Below$WS)



####Boxplot Stuff ####

group_by(PBG_2021_Weight_Below, "WS")

library(ggplot2)

ggplot(PBG_2021_Weight_Below, aes(x = WS, y = `Weight (mg)`, fill = Treatment)) +
  geom_boxplot() +
  labs(title = "Weight by Treatment",
       x = "Treatment",
       y = "Weight (mg)") +
  scale_fill_manual(values = c("blue", "red")) +
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    panel.background = element_rect(fill = "white")
  )

ggsave("myplot.png", dpi = 300)


#### Bargraph ####

ggplot(PBG_2021_Weight_Below, aes(x = Treatment, y = `Weight (mg)`, fill = Treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Weight by Treatment",
       x = "Treatment",
       y = "Weight (mg)") +
  scale_fill_manual(values = c("blue", "red"))

#Ask Kim if these makes sense?

table(PBG_2021_Weight_Below$WS)

# Clean up ####
BurnInfo <- read_excel("YearsSenseBurned.xlsx")


Joined <- full_join(BurnInfo, PBG_2021_Weight_Below)



Joined <- Joined %>% 
  unite("Sample", c("WS", "Trans", "Dist."), sep = "_") %>%
  mutate(Sample = str_remove(Sample, "m")) %>% 
  rename_at('Weight (mg)', ~'Weight')


Joined$Block <- ifelse(grepl("S", Joined $Sample), "North", "South")

Joined <- Joined %>%  unite("Treatment", c("Treatment","SB"), sep = "_")

###Graphs New####

ggplot(Joined, aes(x = Treatment, y = Weight, fill = Treatment)) +
  geom_boxplot() +
  scale_fill_manual(values = c("ABG_0" = "blue", "PBG_0" = "red", "PBG_1" = "red", "PBG_2" = "red")) +
  theme_bw() +
  guides(fill = FALSE) +
  ylab("Weight (mg)")

ggplot(PBG_2021_Weight_Below, aes(x = Treatment, y = `Weight (mg)`, fill = Treatment)) +
  geom_boxplot() +
  scale_fill_manual(values = c("ABG" = "blue", "PBG" = "red")) +
  guides(fill = FALSE) +
  theme(panel.background = element_blank()) +
  theme_bw()


#Stats ####

install.packages("emmeans")

library(nlme) #library for mixed models
library(emmeans) #library for means comparisions after running the models

richModel <- #stores the model output into a named list
  lme(Weight ~ Treatment, #relates richness (or any other dependent variable) to treatment (e.g., ABG and PBG or ABG, PBG0, PBG1, PBG2 for years since burnign)
      data = Joined, #dataset you are analyzing, this must contain all the data (both treatments and all plots)
      random = ~1|Block) #this would be where you'd say north or south unit (which should be a variable in the dataframe)
anova.lme(richModel, type='sequential') #this gives you the ANOVA output from the model, where "sequential" tells it to do a type III anova
emmeans(richModel, pairwise~Treatment, adjust="tukey") #this gives you contrasts (means and confidence intervals) for each possible pairwise comparison of treatments to know whether they are different or not (overlapping confidence intervals means not different)


#### CV ####


# Calculate the CV of weights for each unique watershed
library(dplyr)

# Calculate the CV of weights for each transect within each ABG watershed

PBG_2021_Weight_Below$Block <- ifelse(grepl("S", PBG_2021_Weight_Below$WS), "North", "South")


ABG_CV <- PBG_2021_Weight_Below %>%
  filter(WS %in% c("C1A", "C1SB")) %>%
  group_by(WS) %>%
  summarize(cv_weight = sd(`Weight (mg)`) / mean(`Weight (mg)`) * 100)

# Calculate the CV of weights for each combined PBG cluster
PBG_CV <- PBG_2021_Weight_Below %>%
  filter(WS %in% c("C3A", "C3B", "C3C", "C3SA", "C3SB")) %>%
  group_by(Block) %>%
  summarize(cv_weight = sd(`Weight (mg)`) / mean(`Weight (mg)`) * 100)
