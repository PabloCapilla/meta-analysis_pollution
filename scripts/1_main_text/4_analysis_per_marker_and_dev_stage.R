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
#' Meta-analysis per redox marker and developmental stage
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
phylo_cor <- readRDS("./results/clean_data/data_phylo_cor_20240513.RDS")

##
## effect size data
df00  <- readRDS("./results/clean_data/clean_analysis_20240513.RDS")
head(df00)

#####

##
##### Model 3.1: biomarkers for embryos #####
##

# data for model
df_embryo <- df00 %>% 
  filter(Developmental.Stage == "Embryo")

# model
model31_marker_embryo <- rma.mv(yi=lnRR, 
                        V=lnRR.sv, 
                        mods = ~ Biomarker.Category - 1,
                        random = list(~1 | References,
                                      ~1 | Species,
                                      ~1 | Species_phylo,
                                      ~1 | Observations), 
                        R = list(Species_phylo = phylo_cor),
                        Rscale = "cor",
                        data = df_embryo, 
                        method = "REML", 
                        control=list(optimizer="optimParallel",ncpus=4),
                        sparse = F)
summary(model31_marker_embryo)
round(r2_ml(model = model31_marker_embryo),  digits = 4)   # % var explained by moderator
     
# save model output
saveRDS(object = model31_marker_embryo, file = "./results/models/model31_marker_embryo.RDS")

#####

##
##### Plot model 3.1 #####
##
biomarker_embryo_plot <- orchard_plot(object = model31_marker_embryo, 
                                      mod = "Biomarker.Category", 
                                      group = "References",
                                      trunk.size = 2.0,
                                      cb = FALSE,
                                      xlab = "lnRR", 
                                      transfm = "none") +
  theme(axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_y_continuous(limits = c(-3.75,3.75))

ggsave(filename = "./plots/Figure 2a.pdf", 
       plot = biomarker_embryo_plot, 
       device = "pdf", 
       height = 75, 
       width = 200, 
       units = "mm")

#####

##
##### Model 3.2: biomarkers for tadpole #####
##

# data for model
df_tadpole <- df00 %>% 
  filter(Developmental.Stage == "Tadpole")

# model
model32_marker_Tadpole <- rma.mv(yi=lnRR, 
                                 V=lnRR.sv, 
                                 mods = ~ Biomarker.Category - 1,
                                 random = list(~1 | References,
                                               ~1 | Species,
                                               ~1 | Species_phylo,
                                               ~1 | Observations), 
                                 R = list(Species_phylo = phylo_cor),
                                 Rscale = "cor",
                                 data = df_tadpole, 
                                 method = "REML", 
                                 control=list(optimizer="optimParallel",ncpus=4))
##
## The rma.mv function produces a warning due to ratio of largest to smallest sampling variance extremely large. 
## We have checked that removing the extreme lnRR values that cause the warning and re-running the model produces very similar results.

summary(model32_marker_Tadpole)
round(r2_ml(model = model32_marker_Tadpole),  digits = 4)   # % var explained by moderator

# save model output
saveRDS(object = model32_marker_Tadpole, file = "./results/models/model32_marker_Tadpole.RDS")

#####

##
##### Plot model 3.2 #####
##
biomarker_tadpole_plot <- orchard_plot(object = model32_marker_Tadpole, 
                                          mod = "Biomarker.Category", 
                                          group = "References",
                                          trunk.size = 2.0,
                                          cb = FALSE,
                                          xlab = "lnRR", 
                                          transfm = "none") +
  theme(axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_y_continuous(limits = c(-3.75,3.75))

ggsave(filename = "./plots/Figure 2b.pdf", 
       plot = biomarker_tadpole_plot, 
       device = "pdf", 
       height = 75, 
       width = 200, 
       units = "mm")

#####

##
##### Model 3.3: biomarkers for adults #####
##

# data for model
df_adult <- df00 %>% 
  filter(Developmental.Stage == "Adult")

# model
model33_marker_adult <- rma.mv(yi=lnRR, 
                                V=lnRR.sv, 
                                mods = ~ Biomarker.Category - 1,
                                random = list(~1 | References,
                                              ~1 | Species,
                                              ~1 | Species_phylo,
                                              ~1 | Observations), 
                                R = list(Species_phylo = phylo_cor),
                                Rscale = "cor",
                                data = df_adult, 
                                method = "REML", 
                                control=list(optimizer="optimParallel",ncpus=4),
                                sparse = F)
##
## The rma.mv function produces a warning due to ratio of largest to smallest sampling variance extremely large. 
## We have checked that removing the extreme lnRR values that cause the warning and re-running the model produces very similar results.

summary(model33_marker_adult)
round(r2_ml(model = model33_marker_adult), digits = 4)   # % var explained by moderator
     
# save model output
saveRDS(object = model33_marker_adult, file = "./results/models/model33_marker_adult.RDS")


##
##### Plot model 3.3 #####
##
biomarker_adult_plot <- orchard_plot(object = model33_marker_adult, 
                                          mod = "Biomarker.Category", 
                                         group = "References",
                                         trunk.size = 2,
                                         cb = FALSE,
                                         xlab = "lnRR", 
                                         transfm = "none") +
  theme(axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_y_continuous(limits = c(-3.75,3.75))

ggsave(filename = "./plots/Figure 2c.pdf", 
       plot = biomarker_adult_plot, 
       device = "pdf", 
       height = 75, 
       width = 200, 
       units = "mm")

#####
