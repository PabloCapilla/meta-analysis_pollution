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
#' Meta-analysis per redox marker, developmental stage and pollutant type (organic vs. inorganic)
#' 
##

##
##### libraries & functions ####
##
pacman::p_load(optimParallel, openxlsx, dplyr, tidyr, metafor, ggplot2, wesanderson, orchaRd, ggpubr) 
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
## data summaries 
table(df00$Developmental.Stage, df00$Pollutant.Class_3)

##
##### Models for Organic compound #####
##

##
##### embryos #####
df_organic_embryos <- df00 %>% 
  filter(Pollutant.Class_3 == "Organic compound") %>% 
  filter(Developmental.Stage == "Embryo")

# model
model_organic_embryos <- rma.mv(yi=lnRR, 
                                V=lnRR.sv, 
                                mods = ~ Biomarker.Category - 1, 
                                random = list(~1 | References,
                                              ~1 | Species,
                                              ~1 | Species_phylo,
                                              ~1 | Observations), 
                                R = list(Species_phylo = phylo_cor),
                                Rscale = "cor",
                                data = df_organic_embryos, 
                                method = "REML", 
                                sparse = F, 
                                verbose = T)
summary(model_organic_embryos)

# plot
plot_organic_embryos <- orchard_plot(object = model_organic_embryos, 
                                     mod = "Biomarker.Category", 
                                     group = "References",
                                     trunk.size = 10,
                                     cb = FALSE,
                                     xlab = "lnRR", 
                                     transfm = "none") +
  theme(axis.title = element_text(size = 10),
        axis.text.y = element_text(angle = 45),
        axis.text.x = element_text(size = 10)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) +
  labs(title = "Embryos") +
  scale_y_continuous(limits = c(-4.0,4.0))

#####

##
##### tadpoles #####
df_organic_tadpoles <- df00 %>% 
  filter(Pollutant.Class_3 == "Organic compound") %>% 
  filter(Developmental.Stage == "Tadpole")

# model
model_organic_tadpoles <- rma.mv(yi=lnRR, 
                                 V=lnRR.sv, 
                                 mods = ~ Biomarker.Category - 1, 
                                 random = list(~1 | References,
                                               ~1 | Species,
                                               ~1 | Species_phylo,
                                               ~1 | Observations), 
                                 R = list(Species_phylo = phylo_cor),
                                 Rscale = "cor",
                                 data = df_organic_tadpoles, 
                                 method = "REML", 
                                 sparse = F, 
                                 verbose = T)
summary(model_organic_tadpoles)

# plot
plot_organic_tadpoles <- orchard_plot(object = model_organic_tadpoles, 
                                      mod = "Biomarker.Category", 
                                      group = "References",
                                      trunk.size = 10,
                                      cb = FALSE,
                                      xlab = "lnRR", 
                                      transfm = "none") +
  theme(axis.title = element_text(size = 10),
        axis.text.y = element_text(angle = 45),
        axis.text.x = element_text(size = 10)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) +
  labs(title = "Tadpoles") +
  scale_y_continuous(limits = c(-4.0,4.0))

#####

##
##### adults #####

# data adults
df_organic_adult <- df00 %>% 
  filter(Pollutant.Class_3 == "Organic compound") %>% 
  filter(Developmental.Stage == "Adult")

# model
model_organic_adult <- rma.mv(yi=lnRR, 
                              V=lnRR.sv, 
                              mods = ~ Biomarker.Category - 1, 
                              random = list(~1 | References,
                                            ~1 | Species,
                                            ~1 | Species_phylo,
                                            ~1 | Observations), 
                              R = list(Species_phylo = phylo_cor),
                              Rscale = "cor",
                              data = df_organic_adult, 
                              method = "REML", 
                              sparse = F, 
                              verbose = T)
summary(model_organic_adult)

# plot
plot_organic_adult <- orchard_plot(object = model_organic_adult, 
                                   mod = "Biomarker.Category", 
                                   group = "References",
                                   trunk.size = 10,
                                   cb = FALSE,
                                   xlab = "lnRR", 
                                   transfm = "none") +
  theme(axis.title = element_text(size = 10),
        axis.text.y = element_text(angle = 45),
        axis.text.x = element_text(size = 10)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) +
  labs(title = "Adults") +
  scale_y_continuous(limits = c(-4.0,4.0))

#####

##
##### Panel for Organic compounds #####
##
organic_panel <- ggarrange(plot_organic_embryos, plot_organic_tadpoles, plot_organic_adult,
                           ncol = 1, nrow = 3)

ggsave(filename = "./plots/Figure 3a.pdf", 
       plot = organic_panel, 
       device = "pdf", 
       height = 180, 
       width = 195, 
       units = "mm")

#####

##
##### Models for Inorganic compound #####
##

##
##### embryos #####
df_inorganic_embryos <- df00 %>% 
  filter(Pollutant.Class_3 == "Inorganic compound") %>% 
  filter(Developmental.Stage == "Embryo")

# model
model_inorganic_embryos <- rma.mv(yi=lnRR, 
                                  V=lnRR.sv, 
                                  mods = ~ Biomarker.Category - 1, 
                                  random = list(~1 | References,
                                                ~1 | Species,
                                                ~1 | Species_phylo,
                                                ~1 | Observations), 
                                  R = list(Species_phylo = phylo_cor),
                                  Rscale = "cor",
                                  data = df_inorganic_embryos, 
                                  method = "REML", 
                                  sparse = F, 
                                  verbose = T)
summary(model_inorganic_embryos)

# plot
plot_inorganic_embryos <- orchard_plot(object = model_inorganic_embryos, 
                                       mod = "Biomarker.Category", 
                                       group = "References",
                                       trunk.size = 10,
                                       cb = FALSE,
                                       xlab = "lnRR", 
                                       transfm = "none") +
  theme(axis.title = element_text(size = 10),
        axis.text.y = element_text(angle = 45),
        axis.text.x = element_text(size = 10)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)[c(1,3)]) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)[c(1,3)]) +
  labs(title = "Embryos") +
  scale_y_continuous(limits = c(-4.0,4.0))

#####

##
##### tadpoles #####
df_inorganic_tadpoles <- df00 %>% 
  filter(Pollutant.Class_3 == "Inorganic compound") %>% 
  filter(Developmental.Stage == "Tadpole")

# model
model_inorganic_tadpoles <- rma.mv(yi=lnRR, 
                                   V=lnRR.sv, 
                                   mods = ~ Biomarker.Category - 1, 
                                   random = list(~1 | References,
                                                 ~1 | Species,
                                                 ~1 | Species_phylo,
                                                 ~1 | Observations), 
                                   R = list(Species_phylo = phylo_cor),
                                   Rscale = "cor",
                                   data = df_inorganic_tadpoles, 
                                   method = "REML", 
                                   sparse = F, 
                                   verbose = T)
summary(model_inorganic_tadpoles)

# plot
plot_inorganic_tadpoles <- orchard_plot(object = model_inorganic_tadpoles, 
                                        mod = "Biomarker.Category", 
                                        group = "References",
                                        trunk.size = 10,
                                        cb = FALSE,
                                        xlab = "lnRR", 
                                        transfm = "none") +
  theme(axis.title = element_text(size = 10),
        axis.text.y = element_text(angle = 45),
        axis.text.x = element_text(size = 10)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) +
  labs(title = "Tadpoles") +
  scale_y_continuous(limits = c(-4.0,4.0))

#####

##
##### adults #####

# data adults
df_inorganic_adult <- df00 %>% 
  filter(Pollutant.Class_3 == "Inorganic compound") %>% 
  filter(Developmental.Stage == "Adult")

# model
model_inorganic_adult <- rma.mv(yi=lnRR, 
                                V=lnRR.sv, 
                                mods = ~ Biomarker.Category - 1, 
                                random = list(~1 | References,
                                              ~1 | Species,
                                              ~1 | Species_phylo,
                                              ~1 | Observations), 
                                R = list(Species_phylo = phylo_cor),
                                Rscale = "cor",
                                data = df_inorganic_adult, 
                                method = "REML", 
                                sparse = F, 
                                verbose = T)
summary(model_inorganic_adult)

# plot
plot_inorganic_adult <- orchard_plot(object = model_inorganic_adult, 
                                     mod = "Biomarker.Category", 
                                     group = "References",
                                     trunk.size = 10,
                                     cb = FALSE,
                                     xlab = "lnRR", 
                                     transfm = "none") +
  theme(axis.title = element_text(size = 10),
        axis.text.y = element_text(angle = 45),
        axis.text.x = element_text(size = 10)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) +
  labs(title = "Adults") +
  scale_y_continuous(limits = c(-4.0,4.0))

#####

##
##
##### Panel for Inorganic compounds #####
##
inorganic_panel <- ggarrange(plot_inorganic_embryos, plot_inorganic_tadpoles, plot_inorganic_adult,
                             ncol = 1, nrow = 3)

ggsave(filename = "./plots/Figure 3b.pdf", 
       plot = inorganic_panel, 
       device = "pdf", 
       height = 180, 
       width = 195, 
       units = "mm")

#####
