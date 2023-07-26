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
#' Meta-analysis per redox marker
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
phylo_cor <- readRDS("./results/clean_data/data_phylo_cor_20230724.RDS")

##
## effect size data
df00  <- readRDS("./results/clean_data/clean_analysis_20230724.RDS")
head(df00)

#####

##
##
##### Model 2: Redox Biomarkers #####
##
##
model2_marker <- rma.mv(yi=lnRR, 
                          V=lnRR.sv, 
                          mods = ~ Biomarker.Category - 1,
                          random = list(~1 | References,
                                        ~1 | Species,
                                        ~1 | Species_phylo,
                                        ~1 | Observations), 
                          R = list(Species_phylo = phylo_cor),
                          Rscale = "cor",
                          data = df00, 
                          method = "REML", 
                          control=list(optimizer="optimParallel",ncpus=4),
                          sparse = F)

##
## The rma.mv function produces a warning due to ratio of largest to smallest sampling variance extremely large. 
## We have checked that removing the extreme lnRR values that cause the warning and re-running the model produces very similar results.

summary(model2_marker) # summary
round(r2_ml(model = model2_marker),  digits = 4)  # % var explained by moderator
     

# save model output
saveRDS(object = model2_marker, file = "./results/models/model2_marker.RDS")

#####

##
## re model back if re-starting the script
#overall_model <- readRDS(file = "./results/models/overall_model1.RDS")

##
##
##### Plot model #####
##
##
biomarker_plot <- orchard_plot(object = model2_marker, 
                               mod = "Biomarker.Category", 
                               group = "References",
                               trunk.size = 10,
                               cb = FALSE,
                               xlab = "lnRR", 
                               transfm = "none") +
  theme(axis.title = element_text(size = 15),
        axis.text.x = element_text(size = 15)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) 


ggsave(filename = "./plots/Figure 1b.pdf", 
       plot = biomarker_plot, 
       device = "pdf", 
       height = 100, 
       width = 200, 
       units = "mm")

#####




