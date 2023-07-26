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
#' Meta-analysis across main tissues
#' 
##

##
##
##### libraries & functions ####
##
pacman::p_load(optimParallel, openxlsx, dplyr, tidyr, metafor, ggplot2, ggpubr, wesanderson, orchaRd) 
source("./scripts/0a_R_library/functions.R")

#####

##
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
##
##### Models across pollutants and developmental stage #####
##

#####

##
##### Models for pesticide #####
##

##
## embryos
df_pesticide_embryos <- df00 %>% 
  filter(Pollutant.Class_2b == "Pesticide") %>% 
  filter(Developmental.Stage == "Embryo")

# model
model_pesticide_embryos <- rma.mv(yi=lnRR, 
                                  V=lnRR.sv, 
                                  mods = ~ Biomarker.Category - 1, 
                                  random = list(~1 | References,
                                                ~1 | Species,
                                                ~1 | Species_phylo,
                                                ~1 | Observations), 
                                  R = list(Species_phylo = phylo_cor),
                                  Rscale = "cor",
                                  data = df_pesticide_embryos, 
                                  method = "REML", 
                                  sparse = F, 
                                  verbose = T)
summary(model_pesticide_embryos)

# plot
plot_pesticide_embryos <- orchard_plot(object = model_pesticide_embryos, 
                                       mod = "Biomarker.Category", 
                                       group = "References",
                                       trunk.size = 10,
                                       cb = FALSE,
                                       xlab = "lnRR", 
                                       transfm = "none") +
  theme(axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) +
  labs(title = "Embryos") +
  scale_y_continuous(limits = c(-3.5,3.5))


##
## tadpoles
df_pesticide_tadpoles <- df00 %>% 
  filter(Pollutant.Class_2b == "Pesticide") %>% 
  filter(Developmental.Stage == "Tadpole")

# model
model_pesticide_tadpoles <- rma.mv(yi=lnRR, 
                                   V=lnRR.sv, 
                                   mods = ~ Biomarker.Category - 1, 
                                   random = list(~1 | References,
                                                 ~1 | Species,
                                                 ~1 | Species_phylo,
                                                 ~1 | Observations), 
                                   R = list(Species_phylo = phylo_cor),
                                   Rscale = "cor",
                                   data = df_pesticide_tadpoles, 
                                   method = "REML", 
                                   sparse = F, 
                                   verbose = T)
summary(model_pesticide_tadpoles)

# plot
plot_pesticide_tadpoles <- orchard_plot(object = model_pesticide_tadpoles, 
                                        mod = "Biomarker.Category", 
                                        group = "References",
                                        trunk.size = 10,
                                        cb = FALSE,
                                        xlab = "lnRR", 
                                        transfm = "none") +
  theme(axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) +
  labs(title = "Tadpoles") +
  scale_y_continuous(limits = c(-3.5,3.5))


##
## adults

# data adults
df_pesticide_adult <- df00 %>% 
  filter(Pollutant.Class_2b == "Pesticide") %>% 
  filter(Developmental.Stage == "Adult")

# model
model_pesticide_adult <- rma.mv(yi=lnRR, 
                                V=lnRR.sv, 
                                mods = ~ Biomarker.Category - 1, 
                                random = list(~1 | References,
                                              ~1 | Species,
                                              ~1 | Species_phylo,
                                              ~1 | Observations), 
                                R = list(Species_phylo = phylo_cor),
                                Rscale = "cor",
                                data = df_pesticide_adult, 
                                method = "REML", 
                                sparse = F, 
                                verbose = T)
summary(model_pesticide_adult)

# plot
plot_pesticide_adults <- orchard_plot(object = model_pesticide_adult, 
                                      mod = "Biomarker.Category", 
                                      group = "References",
                                      trunk.size = 10,
                                      cb = FALSE,
                                      xlab = "lnRR", 
                                      transfm = "none") +
  theme(axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) +
  labs(title = "Adults") +
  scale_y_continuous(limits = c(-3.5,3.5))

##
##
##### Panel for pesticides #####
##
##
pesticide_panel <- ggarrange(plot_pesticide_embryos, 
                             plot_pesticide_tadpoles,
                             plot_pesticide_adults, 
                             ncol = 1, nrow = 3)

ggsave(filename = "./plots/Figure S5a - pesticide.pdf", 
       plot = pesticide_panel, 
       device = "pdf", 
       height = 300, 
       width = 175, 
       units = "mm")


#####

#####

##
##### Models for herbicide #####
##


##
## embryos
df_herbicide_embryos <- df00 %>% 
  filter(Pollutant.Class_2b == "Herbicide") %>% 
  filter(Developmental.Stage == "Embryo")

## only 4 observations

##
## tadpoles
df_herbicide_tadpoles <- df00 %>% 
  filter(Pollutant.Class_2b == "Herbicide") %>% 
  filter(Developmental.Stage == "Tadpole")

# model
model_herbicide_tadpoles <- rma.mv(yi=lnRR, 
                                   V=lnRR.sv, 
                                   mods = ~ Biomarker.Category - 1, 
                                   random = list(~1 | References,
                                                 ~1 | Species,
                                                 ~1 | Species_phylo,
                                                 ~1 | Observations), 
                                   R = list(Species_phylo = phylo_cor),
                                   Rscale = "cor",
                                   data = df_herbicide_tadpoles, 
                                   method = "REML", 
                                   sparse = F, 
                                   verbose = T)

summary(model_herbicide_tadpoles)

# plot
plot_herbicide_tadpoles <- orchard_plot(object = model_herbicide_tadpoles, 
                                        mod = "Biomarker.Category", 
                                        group = "References",
                                        trunk.size = 10,
                                        cb = FALSE,
                                        xlab = "lnRR", 
                                        transfm = "none") +
  theme(axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) +
  labs(title = "Tadpoles") +
  scale_y_continuous(limits = c(-3.5,3.5))

##
## adults

# data adults
df_herbicide_adult <- df00 %>% 
  filter(Pollutant.Class_2b == "Herbicide") %>% 
  filter(Developmental.Stage == "Adult")


# model
model_herbicide_adult <- rma.mv(yi=lnRR, 
                                V=lnRR.sv, 
                                mods = ~ Biomarker.Category - 1, 
                                random = list(~1 | References,
                                              ~1 | Species,
                                              ~1 | Species_phylo,
                                              ~1 | Observations), 
                                R = list(Species_phylo = phylo_cor),
                                Rscale = "cor",
                                data = df_herbicide_adult, 
                                method = "REML", 
                                sparse = F, 
                                verbose = T)

summary(model_herbicide_adult)

# plot
plot_herbicide_adult <- orchard_plot(object = model_herbicide_adult, 
                                         mod = "Biomarker.Category", 
                                         group = "References",
                                         trunk.size = 10,
                                         cb = FALSE,
                                         xlab = "lnRR", 
                                         transfm = "none") +
  theme(axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) +
  labs(title = "Adults") +
  scale_y_continuous(limits = c(-3.5,3.5))


##
##
##### Panel for herbicide #####
##
##
herbicide_panel <- ggarrange(plot_herbicide_tadpoles, plot_herbicide_adult,
                             ncol = 1, nrow = 2)

ggsave(filename = "./plots/Figure S5b - herbicide", 
       plot = herbicide_panel, 
       device = "pdf", 
       height = 180, 
       width = 195, 
       units = "mm")

#####

##
##### Models for Metallic elements #####
##

##
## embryos
df_metallic_embryos <- df00 %>% 
  filter(Pollutant.Class_2b == "Metallic elements") %>% 
  filter(Developmental.Stage == "Embryo")

# only 16 estimates from 2 studies

##
## tadpoles
df_metallic_tadpoles <- df00 %>% 
  filter(Pollutant.Class_2b == "Metallic elements") %>% 
  filter(Developmental.Stage == "Tadpole")

# model
model_metallic_tadpoles <- rma.mv(yi=lnRR, 
                                  V=lnRR.sv, 
                                  mods = ~ Biomarker.Category - 1, 
                                  random = list(~1 | References,
                                                ~1 | Species,
                                                ~1 | Species_phylo,
                                                ~1 | Observations), 
                                  R = list(Species_phylo = phylo_cor),
                                  Rscale = "cor",
                                  data = df_metallic_tadpoles, 
                                  method = "REML", 
                                  sparse = F, 
                                  verbose = T)
summary(model_metallic_tadpoles)

# plot
plot_metallic_tadpoles <- orchard_plot(object = model_metallic_tadpoles, 
                                       mod = "Biomarker.Category", 
                                       group = "References",
                                       trunk.size = 10,
                                       cb = FALSE,
                                       xlab = "lnRR", 
                                       transfm = "none") +
  theme(axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) +
  labs(title = "Tadpoles")

##
## adults

# data adults
df_metallic_adult <- df00 %>% 
  filter(Pollutant.Class_2b == "Metallic elements") %>% 
  filter(Developmental.Stage == "Adult")

# model
model_metallic_adult <- rma.mv(yi=lnRR, 
                                V=lnRR.sv, 
                                mods = ~ Biomarker.Category - 1, 
                                random = list(~1 | References,
                                              ~1 | Species,
                                              ~1 | Species_phylo,
                                              ~1 | Observations), 
                                R = list(Species_phylo = phylo_cor),
                                Rscale = "cor",
                                data = df_metallic_adult, 
                                method = "REML", 
                                sparse = F, 
                                verbose = T)
summary(model_metallic_adult)

# plot
plot_metallic_adult <- orchard_plot(object = model_metallic_adult, 
                                         mod = "Biomarker.Category", 
                                        group = "References",
                                        trunk.size = 10,
                                        cb = FALSE,
                                        xlab = "lnRR", 
                                        transfm = "none") +
  theme(axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) +
  labs(title = "Adults")

##
##
##### Panel for Metallic elements #####
##
##
metallic_panel <- ggarrange(plot_metallic_adult, plot_metallic_tadpoles, 
                             ncol = 1, nrow = 2)
ggsave(filename = "./plots/Figure S5c - metallic elements.pdf", 
       plot = metallic_panel, 
       device = "pdf", 
       height = 150, 
       width = 175, 
       units = "mm")
#####

##
##### Models for Inorganic elements #####
##

##
## embryos
df_inorganic_embryos <- df00 %>% 
  filter(Pollutant.Class_2b == "Others_Inorganic_compound") %>% 
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
        axis.text.x = element_text(size = 10)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) +
  labs(title = "Embryos")

##
## tadpoles
df_inorganic_tadpoles <- df00 %>% 
  filter(Pollutant.Class_2b == "Others_Inorganic_compound") %>% 
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
        axis.text.x = element_text(size = 10)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) +
  labs(title = "Tadpoles")

##
## adults

# data adults
df_inorganic_adult <- df00 %>% 
  filter(Pollutant.Class_2b == "Others_Inorganic_compound") %>% 
  filter(Developmental.Stage == "Adult")

# only 9 observations from only 1 study



##
##
##### Panel for Inorganic elements #####
##
##
inorganic_panel <- ggarrange(plot_inorganic_tadpoles, plot_inorganic_embryos, 
                            ncol = 1, nrow = 2)

ggsave(filename = "./plots/Figure S5d - other inorganic elements.pdf.pdf", 
       plot = inorganic_panel, 
       device = "pdf", 
       height = 150, 
       width = 175, 
       units = "mm")
#####
