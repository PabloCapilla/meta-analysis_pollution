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
##### libraries & functions ####
##
pacman::p_load(optimParallel, openxlsx, dplyr, tidyr, metafor, ggplot2, wesanderson, orchaRd, ggpubr) 
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
##### Models for Organic compound #####
##

##
##### embryos #####
table(df00$Pollutant.Class_3)
table(df00$Developmental.Stage)


# model
model_pollutant_dev_state <- rma.mv(yi=lnRR, 
                                V=lnRR.sv, 
                                mods = ~ 
                                  Pollutant.Class_3 +
                                  Developmental.Stage - 1, 
                                random = list(~1 | References,
                                              ~1 | Species,
                                              ~1 | Species_phylo,
                                              ~1 | Observations), 
                                R = list(Species_phylo = phylo_cor),
                                Rscale = "cor",
                                data = df00, 
                                control=list(optimizer="optimParallel",ncpus=8),
                                method = "REML", 
                                sparse = F, 
                                verbose = F)
summary(model_pollutant_dev_state)

# save model output
saveRDS(object = model_pollutant_dev_state, file = "./results/models/s2_model_pollutant_dev_state.RDS")

# plot


HetModel <- orchaRd::mod_results(model_pollutant_dev_state, 
                                 group = "References", 
                                 mod = "Pollutant.Class_3", 
                                 at = list(Developmental.Stage = unique(df00$Developmental.Stage)))

HetModel <- orchaRd::mod_results(model_pollutant_dev_state, 
                                 group = "References", 
                                 mod = "Developmental.Stage", 
                                 at = list(Pollutant.Class_3 = unique(df00$Pollutant.Class_3)))

orchaRd::orchard_plot(HetModel, 
                      xlab = "lnRR",
                      angle = 45, 
                      trunk.size = 2,
                      condition.lab = "Temperature Difference") + 
  theme(legend.direction = "vertical")


plot_pollutant_dev_state <- orchard_plot(object = model_pollutant_dev_state, 
                                     mod = "Biomarker.Category", 
                                     group = "References",
                                     trunk.size = 2,
                                     cb = FALSE,
                                     xlab = "lnRR", 
                                     transfm = "none") +
  theme(axis.title = element_text(size = 10),
        axis.text.y = element_text(angle = 45),
        axis.text.x = element_text(size = 10)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) +
  labs(title = "Embryos") +
  scale_y_continuous(limits = c(-4.0,4.0))

#####
