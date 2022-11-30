###
###
#' 
#' Script for:
#' TITLE
#' Martin et al. 
#' Preprint: 
#' 
#' Latest update: 2022/11/30
#' 
###
###

# Clear memory to make sure there are not files loaded that could cause problems
rm(list=ls())

##
##
##### Script aim: #####
#' This script prepares the original dataset for meta-analysis.
#' 

#####

##
##
##### libraries #####
##
##
pacman::p_load(openxlsx, dplyr, tidyr, 
               metafor, ggplot2, extrafont, 
               RColorBrewer, orchaRd, optimParallel) 
loadfonts()

##
## help functions
source("./scripts/0a_R_library/functions.R")
source("./scripts/0a_R_library/orchard_plot_PCL.R")
source("./scripts/0a_R_library/orchard_plot_PCL_noApples.R")

#####

##
##
##### raw data set #####
##
##

##
## effect size data
df00  <- read.xlsx("./data/Martin_etal_raw_dataset.xlsx",
                   colNames=T,
                   sheet = 1)
head(df00)

## include observation ID
df00$Observations <- 1:nrow(df00)

#####

##
##
##### Calculating adjusted sample sizes #####
##
##

##
## calculate number of rows with the same control data
df00 <- df00 %>% 
  group_by(References, Species, Control_M, Control_SD) %>% 
  mutate(k_shared_control = n()) %>% 
  mutate(Control_N_adj = Control_N / k_shared_control)

# obs with shared controls
df00 %>% 
  filter(k_shared_control > 1) %>% 
  nrow()

# studies with shared controls
df00 %>% 
  filter(k_shared_control > 1) %>% 
  group_by(References) %>% 
  filter(row_number() == 1) %>% 
  nrow()

## double check this calculation
df00 %>% 
  select(References, 
         Species, 
         Control_M, 
         Control_SD, 
         Control_N,
         k_shared_control,
         Control_N_adj) %>% 
  print(n = 20)

## how many observations n control < 1?
df00 %>%
  ungroup() %>% 
  filter(Control_N_adj < 1) %>% 
  summarise(n_obs = n())

# data for effect size calculation
df01 <- df00 %>%
 filter(Control_N_adj >= 1) 

#####

##
##
##### initial sample sizes #####
##
##
nrow(df01) # n observations
length(unique(df01$References)) # n studies
length(unique(df01$Species)) # n species

#####

##
##
##### data cleaning #####
##
##

##
## values with SD = O have been checked by PB ##
##
df01[df01$Treatment_SD == 0, c("References", "Treatment_M", "Treatment_SD", "Treatment_N")]
print(unique(df01[df01$Treatment_SD == 0, c("References")]), n = 50)

df01[df01$Control_SD == 0, c("References", "Control_M", "Control_SD", "Control_N")]
print(unique(df01[df01$Control_SD == 0, c("References")]), n = 50)

##
## remove cases with SD = 0
df01 <- df01 %>% 
  filter(Treatment_SD != 0) %>% 
  filter(Control_SD != 0)

df01[df01$Treatment_SD == 0,]
df01[df01$Control_SD == 0,]

##
## n = 0?
df01[df01$Treatment_N == 0, c("References", "Treatment_M", "Treatment_SD", "Treatment_N")]
df01[df01$Control_N == 0, c("References", "Control_M", "Control_SD", "Control_N")]

#####

##
##
##### Effect size calculation #####
##
##

## SMD
df02 <- as.data.frame(escalc(measure = "SMD", 
                             n2i = Control_N_adj, 
                             sd2i = Control_SD, 
                             m2i = Control_M, 
                             n1i = Treatment_N, 
                             sd1i = Treatment_SD, 
                             m1i = Treatment_M, 
                             var.names = c("SMD", "vSMD"), 
                             data = df01,
                             append=TRUE))
df02 <- df02 %>% 
  filter(!is.na(SMD)) # remove obs with NA SMD

## lnRR
df03 <- as.data.frame(escalc(measure = "ROM", 
                             n2i = Control_N_adj, 
                             sd2i = Control_SD, 
                             m2i = Control_M, 
                             n1i = Treatment_N, 
                             sd1i = Treatment_SD, 
                             m1i = Treatment_M, 
                             var.names = c("lnRR","lnRR.sv"), 
                             data = df02,
                             append=TRUE))

## lnCVR
df04 <- as.data.frame(escalc(measure = "CVR", 
                             n2i = Control_N_adj, 
                             sd2i = Control_SD, 
                             m2i = Control_M, 
                             n1i = Treatment_N, 
                             sd1i = Treatment_SD, 
                             m1i = Treatment_M, 
                             var.names = c("lnCVR","lnCVR.sv"), 
                             data = df03,
                             append=TRUE))


nrow(df04) # n observations
length(unique(df04$References)) # n studies
length(unique(df04$Species)) # n species

##
## dataset with standarsdised effect sizes
df05 <- df04

#####

##
##
##### Removing correlative studies (leaving in Marques et al 2013) #####
##
##
table(df05$Experiment_Venue)

# number of correlative studies (indeed 5)
df05 %>% 
  filter(Experiment_Venue == "Field") %>% 
  group_by(References) %>% 
  filter(row_number() == 1) %>% 
  nrow()

# number of observations from correlative studies
df05 %>% 
  filter(Experiment_Venue == "Field") %>% 
  nrow()

# final dataset
df06 <- df05 %>% 
  filter(Experiment_Venue != "Field")

##
##
##### Summary of sample sizes #####
##
##
nrow(df06) # total of observations
length(unique(df06$References)) # n of studies
length(unique(df06$Species)) # n of species
head(df06)
table(df06$Biomarker_Category)

##
##
##### Observations per groups / levels of moderators #####
##
##

# biomarker / dev stage
table(df06$Biomarker_Category, df06$Developmental_Stage)

df06 %>% 
  group_by(References, Biomarker_Category, Developmental_Stage) %>% 
  filter(row_number() == 1) %>% 
  group_by(Biomarker_Category, Developmental_Stage) %>% 
  summarise(n_obs = n())

# biomarker / dev stage / pollutant
table(df06$Biomarker_Category, df06$Developmental_Stage, df06$Pollutant_Class)

df06 %>% 
  group_by(References, 
           Biomarker_Category, 
           Developmental_Stage, 
           Pollutant_Class) %>% 
  filter(row_number() == 1) %>% 
  group_by(Pollutant_Class, Developmental_Stage, Biomarker_Category) %>% 
  summarise(n_obs = n()) %>% 
  print(n=36)

## testing numbers obtained above
test <- df06 %>% 
  filter(Pollutant_Class == "Metallic elements") %>% 
  filter(Developmental_Stage == "Tadpole") %>% 
  filter(Biomarker_Category == "Enzymatic")
nrow(test)
length(unique(test$References))

#####

##
##
##### Changing species names to match phylo covar matrix and cross-checking #####
##
##
df06[df06$Species == "Lithobates pipiens", "Species"] <- "Rana pipiens" 
df06[df06$Species == "Lithobates catesbeianus", "Species"] <- "Rana catesbeiana" 

##
##
##### Saving data for analysis #####
##
##
saveRDS(object = df06, file = "./data/Martin_etal_analysis_dataset.RDS")





