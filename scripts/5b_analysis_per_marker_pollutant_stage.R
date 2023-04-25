##
##
##### Meta-analysis of Colette's dataset
##### Date: Feb 2022
##### PCL
##
##

rm(list = ls())

##
## Script 5b: meta-analysis per pollutant and developmental stage

##
##### libraries & functions ####
##
pacman::p_load(openxlsx, dplyr, tidyr, 
               metafor, ggplot2, ggpubr, extrafont, 
               RColorBrewer, orchaRd, optimParallel, wesanderson) 
loadfonts()

source("./scripts/0a_R_library/functions.R")
source("./scripts/0a_R_library/orchard_plot_PCL.R")
source("./scripts/0a_R_library/orchard_plot_PCL_noApples.R")

##
##### data #####
##

## 
## phylogeny corr matrix
phylo_cor <- readRDS("./results/clean_data/20221011_data_phylo_cor.RDS")

##
## effect size data
df00  <- readRDS("./results/clean_data/clean_analysis_20221011.RDS")
head(df00)

##
##
##### Changing species names to match phylo covar matrix and cross-checking #####
##
##
df00[df00$Species == "Lithobates pipiens", "Species"] <- "Rana pipiens" 
df00[df00$Species == "Lithobates catesbeianus", "Species"] <- "Rana catesbeiana" 
df00$Species_phylo <- df00$Species # replicating species column for phylo analysis
head(df00)

## matching names between dataset and plylo cor
table(colnames(phylo_cor) %in% df00$Species)
table(colnames(phylo_cor) %in% df00$Species_phylo)

table(df00$Species %in% colnames(phylo_cor))
table(df00$Species_phylo %in% colnames(phylo_cor))


##
##### Models for herbicide #####
##
# data for model
table(df00$Developmental.Stage, df00$Pollutant.Class_2b)

##
## adults

# data adults
df_herbicide_adult <- df00 %>% 
  filter(Pollutant.Class_2b == "Herbicide") %>% 
  filter(Developmental.Stage == "Adult")

table(df_herbicide_adult$Biomarker.Category)

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
                                         data = df_herbicide_adult,
                                         xlab = "lnRR", 
                                         transfm = "none") +
  theme(axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) +
  labs(title = "Adults")

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
                                        data = df_herbicide_tadpoles,
                                        xlab = "lnRR", 
                                        transfm = "none") +
  theme(axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) +
  labs(title = "Tadpoles")


##
## embryos
df_herbicide_embryos <- df00 %>% 
  filter(Pollutant.Class_2b == "Herbicide") %>% 
  filter(Developmental.Stage == "Embryo")

##
##
##### Panel for herbicide #####
##
##
herbicide_panel <- ggarrange(plot_herbicide_adult, plot_herbicide_tadpoles, 
                             ncol = 1, nrow = 2)
ggsave(filename = "./plots/herbicide_dev_panel.pdf", 
       plot = herbicide_panel, 
       device = "pdf", 
       height = 150, 
       width = 175, 
       units = "mm")

#####

##
##### Models for Metallic elements #####
##
# data for model
table(df00$Developmental.Stage, df00$Pollutant.Class_2b)

##
## adults

# data adults
df_metallic_adult <- df00 %>% 
  filter(Pollutant.Class_2b == "Metallic elements") %>% 
  filter(Developmental.Stage == "Adult")

table(df_metallic_adult$Biomarker.Category)

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
                                        data = df_metallic_adult,
                                        xlab = "lnRR", 
                                        transfm = "none") +
  theme(axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) +
  labs(title = "Adults")

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
                                           data = df_metallic_tadpoles,
                                           xlab = "lnRR", 
                                           transfm = "none") +
  theme(axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) +
  labs(title = "Tadpoles")


##
## embryos
df_metallic_embryos <- df00 %>% 
  filter(Pollutant.Class_2b == "Metallic elements") %>% 
  filter(Developmental.Stage == "Embryo")


##
##
##### Panel for Metallic elements #####
##
##
metallic_panel <- ggarrange(plot_metallic_adult, plot_metallic_tadpoles, 
                             ncol = 1, nrow = 2)
ggsave(filename = "./plots/metallic_dev_panel.pdf", 
       plot = metallic_panel, 
       device = "pdf", 
       height = 150, 
       width = 175, 
       units = "mm")
#####

##
##### Models for Inorganic elements #####
##
# data for model
table(df00$Developmental.Stage, df00$Pollutant.Class_2b)

##
## adults

# data adults
df_inorganic_adult <- df00 %>% 
  filter(Pollutant.Class_2b == "Others_Inorganic_compound") %>% 
  filter(Developmental.Stage == "Adult")
table(df_inorganic_adult$Biomarker.Category)

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
                                            data = df_inorganic_tadpoles,
                                            xlab = "lnRR", 
                                            transfm = "none") +
  theme(axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) +
  labs(title = "Tadpoles")


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
                                           data = df_inorganic_embryos,
                                           xlab = "lnRR", 
                                           transfm = "none") +
  theme(axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) +
  labs(title = "Embryos")



##
##
##### Panel for Inorganic elements #####
##
##
inorganic_panel <- ggarrange(plot_inorganic_tadpoles, plot_inorganic_embryos, 
                            ncol = 1, nrow = 2)
ggsave(filename = "./plots/inorganic_dev_panel.pdf", 
       plot = inorganic_panel, 
       device = "pdf", 
       height = 150, 
       width = 175, 
       units = "mm")
#####

##
##### Models for Organic elements #####
##
# data for model
table(df00$Developmental.Stage, df00$Pollutant.Class_2b)

##
## adults

# data adults
df_organic_adult <- df00 %>% 
  filter(Pollutant.Class_2b == "Others_organic_compound") %>% 
  filter(Developmental.Stage == "Adult")
table(df_organic_adult$Biomarker.Category)

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
plot_organic_adults <- orchard_plot(object = model_organic_adult, 
                                            mod = "Biomarker.Category", 
                                        group = "References",
                                        trunk.size = 10,
                                        cb = FALSE,
                                        data = df_organic_adult,
                                        xlab = "lnRR", 
                                        transfm = "none") +
  theme(axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) +
  labs(title = "Adults")


##
## tadpoles
df_organic_tadpoles <- df00 %>% 
  filter(Pollutant.Class_2b == "Others_organic_compound") %>% 
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
                                          data = df_organic_tadpoles,
                                          xlab = "lnRR", 
                                          transfm = "none") +
  theme(axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) +
  labs(title = "Tadpoles")


##
## embryos
df_organic_embryos <- df00 %>% 
  filter(Pollutant.Class_2b == "Others_organic_compound") %>% 
  filter(Developmental.Stage == "Embryo")


##
##
##### Panel for Organic elements #####
##
##
organic_panel <- ggarrange(plot_organic_adults, plot_organic_tadpoles, 
                             ncol = 1, nrow = 2)
ggsave(filename = "./plots/organic_dev_panel.pdf", 
       plot = organic_panel, 
       device = "pdf", 
       height = 150, 
       width = 175, 
       units = "mm")
#####

##
##### Models for Organic elements #####
##
# data for model
table(df00$Developmental.Stage, df00$Pollutant.Class_2b)

##
## adults

# data adults
df_pesticide_adult <- df00 %>% 
  filter(Pollutant.Class_2b == "Pesticide") %>% 
  filter(Developmental.Stage == "Adult")
table(df_pesticide_adult$Biomarker.Category)

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
                                      data = df_pesticide_adult,
                                      xlab = "lnRR", 
                                      transfm = "none") +
  theme(axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) +
  labs(title = "Adults")


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
                                          data = df_pesticide_tadpoles,
                                          xlab = "lnRR", 
                                          transfm = "none") +
  theme(axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) +
  labs(title = "Tadpoles")


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
                                           data = df_pesticide_embryos,
                                           xlab = "lnRR", 
                                           transfm = "none") +
  theme(axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) +
  labs(title = "Embryos")

##
##
##### Panel for pesticides #####
##
##
pesticide_panel <- ggarrange(plot_pesticide_adults, 
                             plot_pesticide_tadpoles, 
                             plot_pesticide_embryos, 
                           ncol = 1, nrow = 3)
ggsave(filename = "./plots/pesticide_dev_panel.pdf", 
       plot = pesticide_panel, 
       device = "pdf", 
       height = 300, 
       width = 175, 
       units = "mm")
#####










