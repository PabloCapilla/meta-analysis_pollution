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
#' Analysis of publication bias
#' 
##

##
##### libraries & functions ####
##
pacman::p_load(openxlsx, dplyr, tidyr, metafor, ggplot2, orchaRd) 
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

## creating new variables to test small-study effects and time lag effects
df00 <- df00 %>% 
  mutate(sqrt_inv_eff_ss = sqrt((1/Treatment.N) + (1/Control.N.adj)),
         publication_year = readr::parse_number(sub(pattern = ' ', 
                                                replacement = '', 
                                                x = substr(df00$References, 
                                                           nchar(df00$References) - 5 + 1, 
                                                           nchar(df00$References)))))
head(df00)
summary(df00$publication_year)

## mean center publication year
df00$publication_year_c <- df00$publication_year - mean(df00$publication_year)
  
#####

##
##### Publication bias model #####
##
##

## small study effects
small_study_model <- rma.mv(yi=lnRR, 
                         V=lnRR.sv, 
                         mods = ~ 1 +
                           sqrt_inv_eff_ss,
                         random = list(~1 | References,
                                       ~1 | Species,
                                       ~1 | Species_phylo,
                                       ~1 | Observations), 
                         R = list(Species_phylo = phylo_cor),
                         Rscale = "cor",
                         data = df00, 
                         method = "REML", 
                         control=list(optimizer="optimParallel",ncpus=6),
                         sparse = F)
##
## The rma.mv function produces a warning due to ratio of largest to smallest sampling variance extremely large. 
## We have checked that removing the extreme lnRR values that cause the warning and re-running the model produces very similar results.

# save model output
saveRDS(object = small_study_model, file = "./results/models/pub_bias_small_study_model.RDS")

summary(small_study_model)

## time lags
time_lag_model <- rma.mv(yi=lnRR, 
                        V=lnRR.sv, 
                        mods = ~ 1 +
                          publication_year_c,
                        random = list(~1 | References,
                                      ~1 | Species,
                                      ~1 | Species_phylo,
                                      ~1 | Observations), 
                        R = list(Species_phylo = phylo_cor),
                        Rscale = "cor",
                        data = df00, 
                        method = "REML", 
                        control=list(optimizer="optimParallel",ncpus=5),
                        sparse = F)
##
## The rma.mv function produces a warning due to ratio of largest to smallest sampling variance extremely large. 
## We have checked that removing the extreme lnRR values that cause the warning and re-running the model produces very similar results.

summary(time_lag_model)

# save model output
saveRDS(object = time_lag_model, file = "./results/models/pub_bias_time_lag_model.RDS")

#####
