##
##
##### Meta-analysis of Colette's dataset
##### Date: Feb 2022
##### PCL
##
##

rm(list = ls())

##
## Script 8: meta-analysis per biomarker (combining enz and non-enz)
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
## phylogeny corr matrix
phylo_cor <- readRDS("./results/clean_data/20220225_data_phylo_cor.RDS")

##
## effect size data
df00  <- readRDS("./results/clean_data/clean_analysis_20220225.RDS")
head(df00)

##
##
##### Changing species names to match phylo covar matrix and cross-checking #####
##
##
df00[df00$Species == "Lithobates pipiens", "Species"] <- "Rana pipiens" 
df00[df00$Species == "Lithobates catesbeianus", "Species"] <- "Rana catesbeiana" 
df00$Species_phylo <- df00$Species # replicating species column for phylo analysis
head(df00)

## matching names between dataset and plylo cor
table(colnames(phylo_cor) %in% df00$Species)
table(colnames(phylo_cor) %in% df00$Species_phylo)

table(df00$Species %in% colnames(phylo_cor))
table(df00$Species_phylo %in% colnames(phylo_cor))


df00$Species_phylo[!df00$Species_phylo %in% colnames(phylo_cor)]

##
##
##### Model 10.0: analysis for main tissues per biomarker #####
##
##

##
## CHECK OBSERVATIONS WITH V LARGE SAMPLING VARIANCE ##
## two observations with very large sampling variance make the model break. 
## I remove then and need to be checked. High SV due to low adj sample size??

df00 %>% 
  filter(lnRR.sv > 10)

## 
## trimmed dataset removing top 10 and bottom 10 obs lnRR.sv observations
df00_trimmed <- df00 %>% 
  filter(lnRR.sv > sort(df00$lnRR.sv)[10]) %>% 
  filter(lnRR.sv < sort(df00$lnRR.sv, decreasing = T)[10])


##
##### Model 10.0.1 - Brain #####
##
table(df00$Biological.Matrix)
df00_brain <- df00 %>% 
  filter(Biological.Matrix == "Brain")

# time tests
system.time(
  model_brain <- rma.mv(yi=lnRR, 
                              V=lnRR.sv, 
                              mods = ~ Biomarker.Category - 1,
                              random = list(~1 | References,
                                            ~1 | Species,
                                            ~1 | Species_phylo,
                                            ~1 | Observations), 
                              R = list(Species_phylo = phylo_cor),
                              Rscale = "cor",
                              data = df00_brain, 
                              method = "REML", 
                              control=list(optimizer="optimParallel",ncpus=3),
                              sparse = F)
)

summary(model_brain) # summary
round(r2_ml(model = model_brain),  # % var explained by moderator
      digits = 4) 


##
##
## Plot brain ##
##
##
library(wesanderson)
names(wes_palettes)
brain_plot <- orchard_plot_PCL(object = model_brain, 
                                   mod = "Biomarker.Category", 
                                   est_point_size = 5,
                                   alpha = 0.5,
                                   cb = FALSE,
                                   xlab = "lnRR", 
                                   ylab = "Biomarker",
                                   transfm = "none", 
                                   angle = 45) +
  theme(axis.title = element_text("Arial", size = 15),
        axis.text.x = element_text("Arial", size = 15)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) +
  labs(title = "Brain")


ggsave(filename = "./plots/model10.0_brain.jpeg", 
       plot = brain_plot, 
       device = "jpeg", 
       height = 100, 
       width = 200, 
       units = "mm")

##
##### Model 10.0.2 - Heart #####
##
table(df00$Biological.Matrix)
df00_heart <- df00 %>% 
  filter(Biological.Matrix == "Heart")


# time tests
system.time(
  model_heart <- rma.mv(yi=lnRR, 
                        V=lnRR.sv, 
                        mods = ~ Biomarker.Category - 1,
                        random = list(#~1 | References,
                                     # ~1 | Species,
                                      #~1 | Species_phylo,
                                      ~1 | Observations), 
                        #R = list(Species_phylo = phylo_cor),
                        #Rscale = "cor",
                        data = df00_heart, 
                        method = "REML", 
                        control=list(optimizer="optimParallel",ncpus=3),
                        sparse = F)
)

summary(model_heart) # summary
round(r2_ml(model = model_heart),  # % var explained by moderator
      digits = 4) 


##
##
## Plot heart ##
##
##
library(wesanderson)
names(wes_palettes)
heart_plot <- orchard_plot_PCL(object = model_heart, 
                               mod = "Biomarker.Category", 
                               est_point_size = 5,
                               alpha = 0.5,
                               cb = FALSE,
                               xlab = "lnRR", 
                               ylab = "Biomarker",
                               transfm = "none", 
                               angle = 45) +
  theme(axis.title = element_text("Arial", size = 15),
        axis.text.x = element_text("Arial", size = 15)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) +
  labs(title = "heart")


ggsave(filename = "./plots/model10.0_heart.jpeg", 
       plot = heart_plot, 
       device = "jpeg", 
       height = 100, 
       width = 200, 
       units = "mm")

##
##### Model 10.0.3 - kidney #####
##
table(df00$Biological.Matrix)
df00_kidney <- df00 %>% 
  filter(Biological.Matrix == "Kidney")


# time tests
system.time(
  model_kidney <- rma.mv(yi=lnRR, 
                        V=lnRR.sv, 
                        mods = ~ Biomarker.Category - 1,
                        random = list(~1 | References,
                           ~1 | Species,
                          ~1 | Species_phylo,
                          ~1 | Observations), 
                        R = list(Species_phylo = phylo_cor),
                        Rscale = "cor",
                        data = df00_kidney, 
                        method = "REML", 
                        control=list(optimizer="optimParallel",ncpus=3),
                        sparse = F)
)

summary(model_kidney) # summary
round(r2_ml(model = model_kidney),  # % var explained by moderator
      digits = 4) 


##
##
## Plot kidney ##
##
##
library(wesanderson)
names(wes_palettes)
kidney_plot <- orchard_plot_PCL(object = model_kidney, 
                               mod = "Biomarker.Category", 
                               est_point_size = 5,
                               alpha = 0.5,
                               cb = FALSE,
                               xlab = "lnRR", 
                               ylab = "Biomarker",
                               transfm = "none", 
                               angle = 45) +
  theme(axis.title = element_text("Arial", size = 15),
        axis.text.x = element_text("Arial", size = 15)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) +
  labs(title = "kidney")


ggsave(filename = "./plots/model10.3_kidney.jpeg", 
       plot = kidney_plot, 
       device = "jpeg", 
       height = 100, 
       width = 200, 
       units = "mm")


##
##### Model 10.0.4 - liver #####
##
table(df00$Biological.Matrix)
df00_liver <- df00 %>% 
  filter(Biological.Matrix == "Liver")


# time tests
system.time(
  model_liver <- rma.mv(yi=lnRR, 
                         V=lnRR.sv, 
                         mods = ~ Biomarker.Category - 1,
                         random = list(~1 | References,
                                       ~1 | Species,
                                       ~1 | Species_phylo,
                                       ~1 | Observations), 
                         R = list(Species_phylo = phylo_cor),
                         Rscale = "cor",
                         data = df00_liver, 
                         method = "REML", 
                         control=list(optimizer="optimParallel",ncpus=3),
                         sparse = F)
)

summary(model_liver) # summary
round(r2_ml(model = model_liver),  # % var explained by moderator
      digits = 4) 


##
##
## Plot liver ##
##
##
library(wesanderson)
names(wes_palettes)
liver_plot <- orchard_plot_PCL(object = model_liver, 
                                mod = "Biomarker.Category", 
                                est_point_size = 5,
                                alpha = 0.5,
                                cb = FALSE,
                                xlab = "lnRR", 
                                ylab = "Biomarker",
                                transfm = "none", 
                                angle = 45) +
  theme(axis.title = element_text("Arial", size = 15),
        axis.text.x = element_text("Arial", size = 15)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) +
  labs(title = "liver")


ggsave(filename = "./plots/model10.3_liver.jpeg", 
       plot = liver_plot, 
       device = "jpeg", 
       height = 100, 
       width = 200, 
       units = "mm")

##
##### Model 10.0.5 - muscle #####
##
table(df00$Biological.Matrix)
df00_muscle <- df00 %>% 
  filter(Biological.Matrix == "Muscle")


# time tests
system.time(
  model_muscle <- rma.mv(yi=lnRR, 
                        V=lnRR.sv, 
                        mods = ~ Biomarker.Category - 1,
                        random = list(~1 | References,
                                      ~1 | Species,
                                      ~1 | Species_phylo,
                                      ~1 | Observations), 
                        R = list(Species_phylo = phylo_cor),
                        Rscale = "cor",
                        data = df00_muscle, 
                        method = "REML", 
                        control=list(optimizer="optimParallel",ncpus=3),
                        sparse = F)
)

summary(model_muscle) # summary
round(r2_ml(model = model_muscle),  # % var explained by moderator
      digits = 4) 


##
##
## Plot muscle ##
##
##
library(wesanderson)
names(wes_palettes)
muscle_plot <- orchard_plot_PCL(object = model_muscle, 
                               mod = "Biomarker.Category", 
                               est_point_size = 5,
                               alpha = 0.5,
                               cb = FALSE,
                               xlab = "lnRR", 
                               ylab = "Biomarker",
                               transfm = "none", 
                               angle = 45) +
  theme(axis.title = element_text("Arial", size = 15),
        axis.text.x = element_text("Arial", size = 15)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) +
  labs(title = "muscle")


ggsave(filename = "./plots/model10.3_muscle.jpeg", 
       plot = muscle_plot, 
       device = "jpeg", 
       height = 100, 
       width = 200, 
       units = "mm")

##
##### Model 10.0.6 - Whole body #####
##
table(df00$Biological.Matrix)
df00_whole <- df00 %>% 
  filter(Biological.Matrix == "Whole Body")


# time tests
system.time(
  model_whole <- rma.mv(yi=lnRR, 
                         V=lnRR.sv, 
                         mods = ~ Biomarker.Category - 1,
                         random = list(~1 | References,
                                       ~1 | Species,
                                       ~1 | Species_phylo,
                                       ~1 | Observations), 
                         R = list(Species_phylo = phylo_cor),
                         Rscale = "cor",
                         data = df00_whole, 
                         method = "REML", 
                         control=list(optimizer="optimParallel",ncpus=3),
                         sparse = F)
)

summary(model_whole) # summary
round(r2_ml(model = model_whole),  # % var explained by moderator
      digits = 4) 


##
##
## Plot whole ##
##
##
library(wesanderson)
names(wes_palettes)
whole_plot <- orchard_plot_PCL(object = model_whole, 
                                mod = "Biomarker.Category", 
                                est_point_size = 5,
                                alpha = 0.5,
                                cb = FALSE,
                                xlab = "lnRR", 
                                ylab = "Biomarker",
                                transfm = "none", 
                                angle = 45) +
  theme(axis.title = element_text("Arial", size = 15),
        axis.text.x = element_text("Arial", size = 15)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) +
  labs(title = "whole")


ggsave(filename = "./plots/model10.3_whole.jpeg", 
       plot = whole_plot, 
       device = "jpeg", 
       height = 100, 
       width = 200, 
       units = "mm")





