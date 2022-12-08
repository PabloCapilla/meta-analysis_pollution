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
#' Test for publication bias
#' 

#####

##
##
##### libraries & functions ####
##
##
#devtools::install_github("daniel1noble/orchaRd", force = TRUE)
library(orchaRd)

pacman::p_load(openxlsx, dplyr, tidyr, wesanderson,
               metafor, ggplot2, extrafont, 
               RColorBrewer, optimParallel) 
loadfonts()

source("./scripts/0a_R_library/functions.R")
source("./scripts/0a_R_library/orchard_plot_PCL.R")
source("./scripts/0a_R_library/orchard_plot_PCL_noApples.R")

#####

##
##
##### dataset #####
##
##

## 
## phylogeny correlation matrix
phylo_cor <- readRDS("./data/phylogenetic_correlation_matrix.RDS")

##
## effect size data
df00  <- readRDS("./data/Martin_etal_analysis_dataset.RDS")
df00$Species_phylo <- df00$Species
head(df00)

#####

mean(as.numeric(df00$year_publication))

##
## data for model
df_lnRR_pb <- df00 %>% 
  filter(!is.na(lnRR)) %>% 
  filter(!is.na(lnRR.sv)) %>% 
  mutate(s_weights = 1/lnRR.sv,
         sqrt_inv_eff_ss = sqrt((1/Treatment_N) + (1/Control_N_adj)),
         inv_eff_ss = (1/Treatment_N) + (1/Control_N_adj),
         publication_year_c = as.numeric(year_publication) - mean(as.numeric(year_publication)))

##
##
##### Publication bias #####
##
##

## small-study effects
model_pub_bias1 <- rma.mv(yi = lnRR, 
                          V = lnRR.sv, 
                          mods = ~ 1 +
                            sqrt_inv_eff_ss, 
                          random = list(~1 | References,
                                        ~1 | Species,
                                        ~1 | Species_phylo,
                                        ~1 | Observations), 
                          R = list(Species_phylo = phylo_cor),
                          data=df_lnRR_pb, 
                          method="REML")
saveRDS(object = model_pub_bias1, file = "./results/model12_pub_bias1.RDS")
model_pub_bias1 <- readRDS("./results/model12_pub_bias1.RDS")          

summary(model_pub_bias1)
r2_ml(model_pub_bias1)


## time lag effects
model_pub_bias2 <- rma.mv(yi = lnRR, 
                          V = lnRR.sv, 
                          mods = ~ 1 +
                            publication_year_c, 
                          random = list(~1 | References,
                                        ~1 | Species,
                                        ~1 | Species_phylo,
                                        ~1 | Observations), 
                          R = list(Species_phylo = phylo_cor),
                          data=df_lnRR_pb, 
                          control=list(optimizer="optimParallel",ncpus=4),
                          sparse = F,
                          method="REML")

saveRDS(object = model_pub_bias2, file = "./results/model13_pub_bias2.RDS")
model_pub_bias2 <- readRDS("./results/model13_pub_bias2.RDS")          

summary(model_pub_bias2)
r2_ml(model_pub_bias2)

