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
#' Overall meta-analysis
#' 

#####

##
##
##### libraries & functions ####
##
##
#devtools::install_github("daniel1noble/orchaRd", force = TRUE)
library(orchaRd)

pacman::p_load(openxlsx, dplyr, tidyr, 
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

##
##
##### Overall model : MODEL 1 #####
##
##

##
## read models back
#overall_model1 <- readRDS("./results/overall_model1.RDS")

##
## running model
overall_model1 <- rma.mv(yi=lnRR, 
                          V=lnRR.sv, 
                          mods = ~ 1,
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

summary(overall_model1) # summary
round(i2_ml(model = overall_model1), digits = 4) # heterogeneity

# save model output
saveRDS(object = overall_model1, file = "./results/overall_model1.RDS")

#####

##
##
##### Overall model : MODEL 1 - removing observations with extreme sv values #####
##
##

##
## read models back
overall_model1b <- readRDS("./results/overall_model1b.RDS") # warning not present

##
## data removing extreme lnRR.sv values (including central 95% of the data distribution)
df00_trimmed <- df00 %>% 
  filter(lnRR.sv > quantile(df00$lnRR.sv, probs = 0.025)) %>% 
  filter(lnRR.sv < quantile(df00$lnRR.sv, probs = 0.975))

overall_model1b <- rma.mv(yi=lnRR, 
                         V=lnRR.sv, 
                         mods = ~ 1,
                         random = list(~1 | References,
                                       ~1 | Species,
                                       ~1 | Species_phylo,
                                       ~1 | Observations), 
                         R = list(Species_phylo = phylo_cor),
                         Rscale = "cor",
                         data = df00_trimmed, 
                         method = "REML", 
                         control=list(optimizer="optimParallel",ncpus=4),
                         sparse = F)

##
## save model output
saveRDS(object = overall_model1b, file = "./results/overall_model1b.RDS")

##
## model results
summary(overall_model1b) # summary
round(i2_ml(model = overall_model1b), digits = 4) # heterogeneity

#####

##
##
##### Plot Overall model #####
##
##
overall_plot <- orchard_plot(object = overall_model1, 
                             data = df00,
                             trunk.size = 10, 
                             branch.size = 1.2, 
                             twig.size = 0.25,
                             mod = "1", 
                             group = "References",
                             alpha = 0.25,
                             xlab = "lnRR", 
                             transfm = "none") +
  theme(axis.title = element_text("Arial", size = 12),
        axis.text.y = element_text("Arial", size = 10),
        axis.text.x = element_text("Arial", size = 10), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks =  seq(-4, 4, 0.5), labels = seq(-4, 4, 0.5)) +
  scale_fill_manual(values = "#bdbdbd") +
  scale_color_manual(values = "#bdbdbd") 

##
## save model
ggsave(filename = "./plots/MODEL1_OVERALL_META-ANALYSIS.png", 
       plot = overall_plot, 
       device = "png", 
       height = 100, 
       width = 150, 
       units = "mm")

#####



