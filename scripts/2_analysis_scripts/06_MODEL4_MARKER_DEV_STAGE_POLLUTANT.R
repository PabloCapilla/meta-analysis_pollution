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
#' Meta-analysis including redox marker as moderator per developmental stage and pollutant (HERBICIDE)
#' 

#####

##
##
##### libraries & functions ####
##
##
pacman::p_load(openxlsx, dplyr, tidyr, ggpubr,
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
##### Models for herbicide #####
##
##

##
## embryos
df_herbicide_embryos <- df00 %>% 
  filter(Pollutant_Class == "Herbicide") %>% 
  filter(Developmental_Stage == "Embryo")
nrow(df_herbicide_embryos)

##
## tadpoles
df_herbicide_tadpoles <- df00 %>% 
  filter(Pollutant_Class == "Herbicide") %>% 
  filter(Developmental_Stage == "Tadpole")

# model
model4.1_herbicide_tadpoles <- rma.mv(yi=lnRR, 
                                   V=lnRR.sv, 
                                   mods = ~ Biomarker_Category - 1, 
                                   random = list(~1 | References,
                                                 ~1 | Species,
                                                 ~1 | Species_phylo,
                                                 ~1 | Observations), 
                                   R = list(Species_phylo = phylo_cor),
                                   Rscale = "cor",
                                   data = df_herbicide_tadpoles, 
                                   method = "REML", 
                                   sparse = F, 
                                   verbose = F)
summary(model4.1_herbicide_tadpoles)


# plot
plot_herbicide_tadpoles <- orchard_plot(object = model4.1_herbicide_tadpoles, 
                                        data = df_herbicide_tadpoles,
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
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) +
  labs(title = "Tadpoles")


##
## adults

# data adults
df_herbicide_adult <- df00 %>% 
  filter(Pollutant_Class == "Herbicide") %>% 
  filter(Developmental_Stage == "Adult")


# model
model4.2_herbicide_adult <- rma.mv(yi=lnRR, 
                                V=lnRR.sv, 
                                mods = ~ Biomarker_Category - 1, 
                                random = list(~1 | References,
                                              ~1 | Species,
                                              ~1 | Species_phylo,
                                              ~1 | Observations), 
                                R = list(Species_phylo = phylo_cor),
                                Rscale = "cor",
                                data = df_herbicide_adult, 
                                method = "REML", 
                                sparse = F, 
                                verbose = F)
summary(model4.2_herbicide_adult)

# plot
plot_herbicide_adult <- orchard_plot(object = model4.2_herbicide_adult, 
                                     data = df_herbicide_adult,
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
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) +
  labs(title = "Adults")



##
##
##### Panel for herbicide #####
##
##
herbicide_panel <- ggarrange(plot_herbicide_adult, plot_herbicide_tadpoles, 
                             ncol = 1, nrow = 2)
ggsave(filename = "./plots/MODEL4_HERBICIDES_PANEL.png", 
       plot = herbicide_panel, 
       device = "png", 
       height = 150, 
       width = 175, 
       units = "mm")

#####

