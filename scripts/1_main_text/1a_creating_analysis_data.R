###
###
#' 
#' Script for:
#' The impact of chemical pollution across major life transitions: a meta-analysis on oxidative stress in amphibians
#' Colette Martin, Pablo Capilla-Lasheras, Pat Monaghan, Pablo Burraco
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
##### libraries & functions ####
##
pacman::p_load(openxlsx, dplyr, tidyr, metafor, ggplot2) 
source("./scripts/0a_R_library/functions.R")

#####

##
##### data #####
##

##
## effect size data
df00  <- read.xlsx("./data/Full_dataset_revision.xlsx",
                   colNames=T,
                   sheet = 1)
head(df00)

##
## remove studies that can't be included in the analysis (see reasons in the data table)
df00 <- df00 %>% 
  filter(Included.in.the.analysis == 'Yes')


#####

##
##### Calculating adjusted sample sizes #####
##
df00 <- df00 %>% 
  group_by(References, Species, Control.M, Control.SD) %>% 
  mutate(k_shared_control = n()) %>% 
  mutate(Control.N.adj = Control.N / k_shared_control)

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

#####

##
##### Observations with SD = O? - yes, these have been checked and they are not mistakes, they are real data #####
##
df01[df01$Treatment.SD == 0, "Experiment.Venue"]
df01[df01$Control.SD == 0, "Experiment.Venue"]
df01[df01$Treatment.SD == 0 | df01$Control.SD == 0, "Experiment.Venue"]

df01[df01$Treatment.SD == 0,"Treatment.SD"] <- 0.01
df01[df01$Control.SD == 0,"Control.SD"] <- 0.01

df01[df01$Treatment.SD == 0,]
df01[df01$Control.SD == 0,]

#####

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
                             data = df02,
                             append=TRUE))
head(df03)

#####

##
##### Removing correlative studies (leaving in Marques et al 2013, which is a field experiment) #####
##
table(df03$Experiment.Venue)

# number of correlative studies
df03 %>% 
  filter(Experiment.Venue == "Field") %>% 
  group_by(References) %>% 
  filter(row_number() == 1) %>% 
  nrow()
  
# number of observations from correlative studies
df03 %>% 
  filter(Experiment.Venue == "Field") %>% 
  nrow()

# final dataset
df04 <- df03 %>% 
  filter(Experiment.Venue != "Field")

#####

##
##### Summary of sample sizes #####
##
nrow(df04) # total of observations
length(unique(df04$References)) # n of studies
length(unique(df04$Species)) # n of species
head(df04)

#####

##
##### Correcting species names before saving data #####
##
df04[df04$Species == "Lithobates pipiens", "Species"] <- "Rana pipiens" 
df04[df04$Species == "Lithobates catesbeianus", "Species"] <- "Rana catesbeiana" 
df04[df04$Species == "Aquarana catesbeianus", "Species"] <- "Rana catesbeiana" 
df04[df04$Species == "Rana nigromaculata", "Species"] <- 'Pelophylax nigromaculatus'
df04[df04$Species == "Bufo raddei", "Species"] <- 'Strauchbufo raddei' 
df04[df04$Species == "Rana saharica", "Species"] <- "Pelophylax saharicus"   
df04[df04$Species == "Rana ridibunda", "Species"] <- "Pelophylax ridibundus"

df04$Species_phylo <- df04$Species # replicating species column for phylo analysis
head(df04)

#####

##
##
##### Saving data for analysis #####
##
##
df04$Observations <- 1:nrow(df04)
saveRDS(object = df04, file = "./results/clean_data/clean_analysis_20240527.RDS")





