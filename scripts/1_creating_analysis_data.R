###
###
#' 
#' Script for:
#' TITLE
#' Martin et al. 
#' Preprint: 
#' 
#' Latest update: 2022/11/29
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
##
##
##


##
##### libraries #####
##






##
## Script 1: creating data set for analysis

##
##### libraries & functions ####
##
pacman::p_load(openxlsx, dplyr, tidyr, 
               metafor, ggplot2, extrafont, 
               RColorBrewer, orchaRd, optimParallel) 
loadfonts()

source("./scripts/0a_R_library/functions.R")
source("./scripts/0a_R_library/orchard_plot_PCL.R")
source("./scripts/0a_R_library/orchard_plot_PCL_noApples.R")

##
##### data #####
##

##
## effect size data
df00  <- read.xlsx("./data/Full_dataset_CM_PB.xlsx",
                   colNames=T,
                   sheet = 2)
head(df00)
df00$Biomarker.Category_2 <- ifelse(df00$Biomarker.Category == "Indicator", "Indicator", "Antioxidant")

## include observation ID
df00$Observations <- 1:nrow(df00)

##
##### Calculating adjusted sample sizes #####
##
df00 <- df00 %>% 
  group_by(References, Species, Control.M, Control.SD) %>% 
  mutate(k_shared_control = n()) %>% 
  mutate(Control.N.adj = Control.N / k_shared_control)

##
## number of obs and studies with share controls per observations

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
         Control.M, 
         Control.SD, 
         Control.N,
         k_shared_control,
         Control.N.adj) %>% 
  print(n = 20)

## how many observations n control < 1?
df00 %>%
  ungroup() %>% 
  filter(Control.N.adj < 1) %>% 
  summarise(n_obs = n())

# data for effect size calculation
df01 <- df00 %>%
  filter(Control.N.adj >= 1) 

nrow(df01) # n observations
length(unique(df01$References)) # n studies
length(unique(df01$Species)) # n species


##
##
###
##### VALUES WITH SV = O - NEED TO BE CHECKED!!!! #####
##### CHECKED BY PB - I RETAIN THESE OBS #####
###
##
##
df01[df01$Treatment.SD == 0, "Experiment.Venue"]
df01[df01$Control.SD == 0, "Experiment.Venue"]
df01[df01$Treatment.SD == 0 & df01$Control.SD == 0, "Experiment.Venue"]


df01[df01$Treatment.SD == 0,"Treatment.SD"] <- 0.01
df01[df01$Control.SD == 0,"Control.SD"] <- 0.01

df01[df01$Treatment.SD == 0,]
df01[df01$Control.SD == 0,]


#dft <- df01 %>% 
#  filter(Treatment.SD != 0) %>% 
#  filter(Control.SD != 0)

##
##
##### Effect size calculation #####
##
##

## SMD
df02 <- as.data.frame(escalc(measure = "SMD", 
                             n2i = Control.N.adj, 
                             sd2i = Control.SD, 
                             m2i = Control.M, 
                             n1i = Treatment.N, 
                             sd1i = Treatment.SD, 
                             m1i = Treatment.M, 
                             var.names = c("SMD", "vSMD"), 
                             data = df01,
                             append=TRUE))
head(df02)

## lnRR
df03 <- as.data.frame(escalc(measure = "ROM", 
                             n2i = Control.N.adj, 
                             sd2i = Control.SD, 
                             m2i = Control.M, 
                             n1i = Treatment.N, 
                             sd1i = Treatment.SD, 
                             m1i = Treatment.M, 
                             var.names = c("lnRR","lnRR.sv"), 
                             data = df01,
                             append=TRUE))

## lnCVR
df04 <- as.data.frame(escalc(measure = "CVR", 
                             n2i = Control.N.adj, 
                             sd2i = Control.SD, 
                             m2i = Control.M, 
                             n1i = Treatment.N, 
                             sd1i = Treatment.SD, 
                             m1i = Treatment.M, 
                             var.names = c("lnCVR","lnCVR.sv"), 
                             data = df03,
                             append=TRUE))


nrow(df04) # n observations
length(unique(df04$References)) # n studies
length(unique(df04$Species)) # n species

df05 <- df04

##
##
##### Removing correlative studies (leaving in Marques et al 2013) #####
##
##
table(df05$Experiment.Venue)

# number of correlative studies (indeed 5)
df05 %>% 
  filter(Experiment.Venue == "Field") %>% 
  group_by(References) %>% 
  filter(row_number() == 1) %>% 
  nrow()
  
# number of observations from correlative studies
df05 %>% 
  filter(Experiment.Venue == "Field") %>% 
  nrow()

# final dataset
df06 <- df05 %>% 
  filter(Experiment.Venue != "Field")

##
##
##### Summary of sample sizes #####
##
##
nrow(df06) # total of observations
length(unique(df06$References)) # n of studies
length(unique(df06$Species)) # n of species
head(df06)

##
##
##### Observations per groups / levels of moderators #####
##
##

# biomarker / dev stage
table(df06$Biomarker.Category, df06$Developmental.Stage)

df06 %>% 
  group_by(References, Biomarker.Category, Developmental.Stage) %>% 
  filter(row_number() == 1) %>% 
  group_by(Biomarker.Category, Developmental.Stage) %>% 
  summarise(n_obs = n())
  
# biomarker / dev stage / pollutant
table(df06$Biomarker.Category, df06$Developmental.Stage, df06$Pollutant.Class_2b)

df06 %>% 
  group_by(References, 
           Biomarker.Category, 
           Developmental.Stage, 
           Pollutant.Class_2b) %>% 
  filter(row_number() == 1) %>% 
  group_by(Pollutant.Class_2b, Developmental.Stage, Biomarker.Category) %>% 
  summarise(n_obs = n()) %>% 
  print(n=36)

## testing numbers obtained above
test <- df06 %>% 
  filter(Pollutant.Class_2b == "Metallic elements") %>% 
  filter(Developmental.Stage == "Tadpole") %>% 
  filter(Biomarker.Category == "Enzymatic")
nrow(test)
length(unique(test$References))

##
##
##### Saving data for analysis #####
##
##
saveRDS(object = df06, file = "./results/clean_data/clean_analysis_20221011.RDS")





