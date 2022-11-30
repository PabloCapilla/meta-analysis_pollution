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
#' Meta-analysis including redox marker as moderator
#' 

#####

##
##
##### libraries & functions ####
##
##
pacman::p_load(openxlsx, dplyr, tidyr, 
               metafor, ggplot2, extrafont, 
               RColorBrewer, orchaRd, wesanderson, 
               optimParallel) 
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
##### Model 2: Biomarkers #####
##
##
##
## read models back
#overall_model2 <- readRDS("./results/overall_model2.RDS")

##
## running model
overall_model2 <- rma.mv(yi=lnRR, 
                         V=lnRR.sv, 
                         mods = ~ Biomarker_Category - 1,
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

summary(overall_model2) # summary
round(i2_ml(model = overall_model2), digits = 4) # heterogeneity

# save model output
saveRDS(object = overall_model2, file = "./results/overall_model2.RDS")

#####

##
##
##### Overall model : MODEL 2 - removing observations with extreme sv values #####
##
##

##
## read models back
overall_model2b <- readRDS("./results/overall_model2b.RDS") # warning not present

##
## data removing extreme lnRR.sv values (including central 95% of the data distribution)
df00_trimmed <- df00 %>% 
  filter(lnRR.sv > quantile(df00$lnRR.sv, probs = 0.025)) %>% 
  filter(lnRR.sv < quantile(df00$lnRR.sv, probs = 0.975))

overall_model2b <- rma.mv(yi=lnRR, 
                          V=lnRR.sv, 
                          mods = ~ Biomarker_Category - 1,
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
saveRDS(object = overall_model2b, file = "./results/overall_model2b.RDS")

##
## model results
summary(overall_model2b) # summary
round(i2_ml(model = overall_model2b), digits = 4) # heterogeneity

#####

##
##
##### Plot Overall model #####
##
##

##
## preparing colours for plot
names(wes_palettes)
biomarker_plot <- orchard_plot_PCL(object = overall_model2, 
                                   mod = "Biomarker_Category", 
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
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) 


ggsave(filename = "./plots/MODEL2_MARKER.png", 
       plot = biomarker_plot, 
       device = "png", 
       height = 100, 
       width = 200, 
       units = "mm")





