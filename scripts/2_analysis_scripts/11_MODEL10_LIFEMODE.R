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
#' Meta-analysis per life mode
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

##
##
##### Model 10: life mode #####
##
##

# model
model10_lifemode <- rma.mv(yi=lnRR, 
                        V=lnRR.sv, 
                        mods = ~ Life_cycle - 1,
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
# save model output
saveRDS(object = model10_lifemode, file = "./results/model10_lifemode.RDS")
model10_lifemode <- readRDS("./results/model10_lifemode.RDS")          

##
## model results
summary(model10_lifemode) # summary
round(r2_ml(model = model10_lifemode),  # % var explained by moderator
      digits = 4) 



##
##
##### Plot model #####
##
##
lifecycle_plot <- orchard_plot(object = model10_lifemode, 
                            data = df00,
                            trunk.size = 10, 
                            branch.size = 1.2, 
                            twig.size = 0.25,
                            mod = "Life_cycle", 
                            group = "References",
                            alpha = 0.25,
                            xlab = "lnRR", 
                            transfm = "none") +
  theme(axis.title = element_text("Arial", size = 12),
        axis.text.y = element_text("Arial", size = 10, angle = 0),
        axis.text.x = element_text("Arial", size = 10), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks =  seq(-4, 4, 0.5), labels = seq(-4, 4, 0.5)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 2)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 2)) 

ggsave(filename = "./plots/MODEL10_LIFECYCLE.png", 
       plot = lifecycle_plot, 
       device = "png", 
       height = 100, 
       width = 200, 
       units = "mm")





