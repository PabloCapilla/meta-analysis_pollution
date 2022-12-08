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
##### Model 3.1: biomarkers for embryos #####
##
##

# data for model
df_embryo <- df00 %>% 
  filter(Developmental_Stage == "Embryo")

# model
model31_marker_embryo <- rma.mv(yi=lnRR, 
                                V=lnRR.sv, 
                                mods = ~ Biomarker_Category - 1,
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
#saveRDS(object = model31_marker_embryo, file = "./results//model3.1_marker_embryo.RDS")
#model31_marker_embryo <- readRDS(file = "./results/models/model31_marker_embryo.RDS")

summary(model31_marker_embryo)
round(r2_ml(model = model31_marker_embryo),  # % var explained by moderator
      digits = 4) 


# plot
biomarker_embryo_plot <- orchard_plot(object = model31_marker_embryo, 
                                      data = df_embryo,
                                      trunk.size = 10, 
                                      branch.size = 1.2, 
                                      twig.size = 0.25,
                                      mod = "Biomarker_Category", 
                                      group = "References",
                                      alpha = 0.25,
                                      xlab = "lnRR", 
                                      transfm = "none") +
  theme(axis.title = element_text("Arial", size = 12),
        axis.text.y = element_text("Arial", size = 10, angle = 0),
        axis.text.x = element_text("Arial", size = 10), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks =  seq(-4, 4, 0.5), labels = seq(-4, 4, 0.5)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) 

  
ggsave(filename = "./plots/MODEL3.1_BIOMARKER_EMBRYO.png", 
       plot = biomarker_embryo_plot, 
       device = "png", 
       height = 75, 
       width = 200, 
       units = "mm")

#####

##
##
##### Model 3.2: biomarkers for tadpoles #####
##
##

# data for model
df_tadpole <- df00 %>% 
  filter(Developmental_Stage == "Tadpole")

# model
model32_marker_Tadpole <- rma.mv(yi=lnRR, 
                                 V=lnRR.sv, 
                                 mods = ~ Biomarker_Category - 1,
                                 random = list(~1 | References,
                                               ~1 | Species,
                                               ~1 | Species_phylo,
                                               ~1 | Observations), 
                                 R = list(Species_phylo = phylo_cor),
                                 Rscale = "cor",
                                 data = df_tadpole, 
                                 method = "REML", 
                                 control=list(optimizer="optimParallel",ncpus=3),
                                 sparse = F)
#saveRDS(object = model32_marker_Tadpole, file = "./results/model3.2_marker_tadpole.RDS")
#model32_marker_Tadpole <- readRDS("./results/model3.2_marker_tadpole.RDS")          # warning present

summary(model32_marker_Tadpole)
round(r2_ml(model = model32_marker_Tadpole), 
      digits = 4) 


# plot
biomarker_tadpole_plot <- orchard_plot(object = model32_marker_Tadpole, 
                                       data = df_tadpole,
                                       trunk.size = 10, 
                                       branch.size = 1.2, 
                                       twig.size = 0.25,
                                       mod = "Biomarker_Category", 
                                       group = "References",
                                       alpha = 0.25,
                                       xlab = "lnRR", 
                                       transfm = "none") +
  theme(axis.title = element_text("Arial", size = 12),
        axis.text.y = element_text("Arial", size = 10, angle = 0),
        axis.text.x = element_text("Arial", size = 10), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks =  seq(-4, 4, 0.5), labels = seq(-4, 4, 0.5)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) 

ggsave(filename = "./plots/MODEL3.2_BIOMARKER_TADPOLE.png", 
       plot = biomarker_tadpole_plot, 
       device = "png", 
       height = 75, 
       width = 200, 
       units = "mm")

#####

##
##
##### Model 3.3: biomarkers for adults #####
##
##

# data for model
df_adult <- df00 %>% 
  filter(Developmental_Stage == "Adult")

# model
model33_marker_adult <- rma.mv(yi=lnRR, 
                               V=lnRR.sv, 
                               mods = ~ Biomarker_Category - 1,
                               random = list(~1 | References,
                                             ~1 | Species,
                                             ~1 | Species_phylo,
                                             ~1 | Observations), 
                               R = list(Species_phylo = phylo_cor),
                               Rscale = "cor",
                               data = df_adult, 
                               method = "REML", 
                               control=list(optimizer="optimParallel",ncpus=3),
                               sparse = F)
#saveRDS(object = model33_marker_adult, file = "./results/model3.3_marker_adult.RDS")
#model33_marker_adult <- readRDS("./results/model3.3_marker_adult.RDS")

summary(model33_marker_adult)
round(r2_ml(model = model33_marker_adult),  # % var explained by moderator
      digits = 4) 


# plot
biomarker_adult_plot <- orchard_plot(object = model33_marker_adult, 
                                     data = df_adult,
                                     trunk.size = 10, 
                                     branch.size = 1.2, 
                                     twig.size = 0.25,
                                     mod = "Biomarker_Category", 
                                     group = "References",
                                     alpha = 0.25,
                                     xlab = "lnRR", 
                                     transfm = "none") +
  theme(axis.title = element_text("Arial", size = 12),
        axis.text.y = element_text("Arial", size = 10, angle = 0),
        axis.text.x = element_text("Arial", size = 10), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks =  seq(-4, 4, 0.5), labels = seq(-4, 4, 0.5)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) 

ggsave(filename = "./plots/MODEL3.3_BIOMARKER_ADULT.png", 
       plot = biomarker_adult_plot, 
       device = "png", 
       height = 75, 
       width = 200, 
       units = "mm")

=