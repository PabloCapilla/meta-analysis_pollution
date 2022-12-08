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
#' Meta-analysis including redox marker as moderator per developmental stage and pollutant (PESTICIDE)
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
df_embryos <- df00 %>% 
  filter(Pollutant_Class == "Pesticide") %>% 
  filter(Developmental_Stage == "Embryo")
nrow(df_embryos)

# model
model5.1_embryos <- rma.mv(yi=lnRR, 
                                      V=lnRR.sv, 
                                      mods = ~ Biomarker_Category - 1, 
                                      random = list(~1 | References,
                                                    ~1 | Species,
                                                    ~1 | Species_phylo,
                                                    ~1 | Observations), 
                                      R = list(Species_phylo = phylo_cor),
                                      Rscale = "cor",
                                      data = df_embryos, 
                                      method = "REML", 
                                      sparse = F, 
                                      verbose = F)
summary(model5.1_embryos)


# plot
plot_embryos <- orchard_plot(object = model5.1_embryos, 
                                        data = df_embryos,
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
  labs(title = "Embryos")

##
## tadpoles
df_tadpoles <- df00 %>% 
  filter(Pollutant_Class == "Pesticide") %>% 
  filter(Developmental_Stage == "Tadpole")

# model
model5.2_tadpoles <- rma.mv(yi=lnRR, 
                                      V=lnRR.sv, 
                                      mods = ~ Biomarker_Category - 1, 
                                      random = list(~1 | References,
                                                    ~1 | Species,
                                                    ~1 | Species_phylo,
                                                    ~1 | Observations), 
                                      R = list(Species_phylo = phylo_cor),
                                      Rscale = "cor",
                                      data = df_tadpoles, 
                                      method = "REML", 
                                      sparse = F, 
                                      verbose = F)
summary(model5.2_tadpoles)


# plot
plot_tadpoles <- orchard_plot(object = model5.2_tadpoles, 
                                        data = df_tadpoles,
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
df_adult <- df00 %>% 
  filter(Pollutant_Class == "Pesticide") %>% 
  filter(Developmental_Stage == "Adult")


# model
model5.3_adult <- rma.mv(yi=lnRR, 
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
                                   sparse = F, 
                                   verbose = F)
summary(model5.3_adult)

# plot
plot_adult <- orchard_plot(object = model5.3_adult, 
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
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) +
  labs(title = "Adults")



##
##
##### Panel for pesticides #####
##
##
pesticide_panel <- ggarrange(plot_adult, plot_tadpoles, plot_embryos,
                             ncol = 1, nrow = 3)
ggsave(filename = "./plots/MODEL5_PESTICIDES_PANEL.png", 
       plot = pesticide_panel, 
       device = "png", 
       height = 150, 
       width = 175, 
       units = "mm")

