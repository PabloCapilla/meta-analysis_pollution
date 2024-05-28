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
##### Script aim: #####
#' Meta-analysis across main tissues
#' 
##

##
##### libraries & functions ####
##
pacman::p_load(optimParallel, openxlsx, dplyr, tidyr, metafor, ggplot2, wesanderson, orchaRd) 
source("./scripts/0a_R_library/functions.R")

#####

##
##### data #####
##

## 
## phylogeny corr matrix
phylo_cor <- readRDS("./results/clean_data/data_phylo_cor_20240527.RDS")

##
## effect size data
df00  <- readRDS("./results/clean_data/clean_analysis_20240527.RDS")
head(df00)

#####

##
##
##### Analysis across main tissues #####
##
##
table(df00$Biological.Matrix)
df_tissue <- df00 %>% 
  filter(Biological.Matrix == "Brain" |
           Biological.Matrix == "Heart" |
           Biological.Matrix == "Kidney"| 
           Biological.Matrix == "Liver" |
           Biological.Matrix == "Muscle" |
           Biological.Matrix == "Whole Body")

##
## model
model_tissue <- rma.mv(yi=lnRR, 
                       V=lnRR.sv, 
                       mods = ~ Biological.Matrix - 1,
                       random = list(~1 | References,
                                     ~1 | Species,
                                     ~1 | Species_phylo,
                                     ~1 | Observations), 
                       R = list(Species_phylo = phylo_cor),
                       Rscale = "cor",
                       data = df_tissue, 
                       method = "REML", 
                       control=list(optimizer="optimParallel",ncpus=8),
                       sparse = F)
##
## The rma.mv function produces a warning due to ratio of largest to smallest sampling variance extremely large. 
## We have checked that removing the extreme lnRR values that cause the warning and re-running the model produces very similar results.


summary(model_tissue) # summary
round(r2_ml(model = model_tissue),  # % var explained by moderator
      digits = 4) 


# save model output
saveRDS(object = model_tissue, file = "./results/models/model_tissue.RDS")

#####

##
##### Plot model #####
##
biomarker_plot <- orchard_plot(object = model_tissue, 
                               mod = "Biological.Matrix", 
                               group = "References",
                               trunk.size = 2.0,
                               cb = FALSE,
                               xlab = "lnRR", 
                               transfm = "none") +
  theme(axis.title = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15, angle = 40)) 

ggsave(filename = "./plots/Figure S4.pdf", 
       plot = biomarker_plot, 
       device = "pdf", 
       height = 100, 
       width = 200, 
       units = "mm")

ggsave(filename = "./plots/Figure S4.png", 
       plot = biomarker_plot, 
       device = "png", 
       height = 100, 
       width = 200, 
       units = "mm")

#####



