##
##
##### Meta-analysis of Colette's dataset
##### Date: Feb 2022
##### PCL
##
##

rm(list = ls())

##
## Script 5: meta-analysis per biomaker and pollutant

##
##### libraries & functions ####
##
pacman::p_load(openxlsx, dplyr, tidyr, 
               metafor, ggplot2, extrafont, 
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
##
##### Sample size per biomarker and pollutant #####
##
##
table(df00$Pollutant.Class_2b, df00$Biomarker.Category)



##
##
##
##### Model 4: Biomarkers and pollutants #####
##
##

##
##### Model 4.1: biomarkers for herbicide #####

# data for model
df_herbicide <- df00 %>% 
  filter(Pollutant.Class_2b == "Herbicide")
table(df_herbicide$Developmental.Stage)


# model
model41_marker_herbicide <- rma.mv(yi=lnRR, 
                                   V=lnRR.sv, 
                                   mods = ~ Biomarker.Category - 1, 
                                   random = list(~1 | References,
                                                 ~1 | Species,
                                                 ~1 | Species_phylo,
                                                 ~1 | Observations), 
                                   R = list(Species_phylo = phylo_cor),
                                   Rscale = "cor",
                                   data = df_herbicide, 
                                   method = "REML", 
                                   #control=list(optimizer="optimParallel",ncpus=2),
                                   sparse = F, 
                                   verbose = T)

#saveRDS(object = model41_marker_herbicide, file = "./results/models/model41_marker_herbicide.RDS")
model41_marker_herbicide <- readRDS(file = "./results/models/model41_marker_herbicide.RDS")

summary(model41_marker_herbicide)
round(r2_ml(model = model41_marker_herbicide),  # % var explained by moderator
      digits = 4) 


# plot
biomarker_herbicide_plot <- orchard_plot(object = model41_marker_herbicide, 
                                          mod = "Biomarker.Category", 
                                          group = "References",
                                          trunk.size = 10,
                                          cb = FALSE,
                                          data = df_herbicide,
                                          xlab = "lnRR", 
                                          transfm = "none") +
  theme(axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) 

ggsave(filename = "./plots/model41_marker_herbicide.pdf", 
       plot = biomarker_herbicide_plot, 
       device = "pdf", 
       height = 75, 
       width = 200, 
       units = "mm")


##
##### Model 4.2: biomarkers for Metallic elements #####

# data for model
df_Metallic_elements <- df00 %>% 
  filter(Pollutant.Class_2b == "Metallic elements")
table(df_Metallic_elements$Developmental.Stage)


# model
model42_marker_Metallic_elements <- rma.mv(yi=lnRR, 
                                   V=lnRR.sv, 
                                   mods = ~ Biomarker.Category - 1, 
                                   random = list(~1 | References,
                                                 ~1 | Species,
                                                 ~1 | Species_phylo,
                                                 ~1 | Observations), 
                                   R = list(Species_phylo = phylo_cor),
                                   Rscale = "cor",
                                   data = df_Metallic_elements, 
                                   method = "REML", 
                                   #control=list(optimizer="optimParallel",ncpus=2),
                                   sparse = F, 
                                   verbose = T)

#saveRDS(object = model42_marker_Metallic_elements, file = "./results/models/model42_marker_Metallic_elements.RDS")
model42_marker_Metallic_elements<- readRDS(file = "./results/models/model42_marker_Metallic_elements.RDS")

summary(model42_marker_Metallic_elements)
round(r2_ml(model = model42_marker_Metallic_elements),  # % var explained by moderator
      digits = 4) 


# plot
marker_Metallic_elements_plot <- orchard_plot(object = model42_marker_Metallic_elements, 
                                             mod = "Biomarker.Category", 
                                             group = "References",
                                             trunk.size = 10,
                                             cb = FALSE,
                                             data = df_Metallic_elements,
                                             xlab = "lnRR", 
                                             transfm = "none") +
  theme(axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) 

ggsave(filename = "./plots/model42_marker_Metallic_elements.pdf", 
       plot = marker_Metallic_elements_plot, 
       device = "pdf", 
       height = 100, 
       width = 200, 
       units = "mm")




##### Model 4.3: biomarkers for Others Inorganic compound #####

# data for model
df_inorganic <- df00 %>% 
  filter(Pollutant.Class_2b == "Others_Inorganic_compound")
table(df_inorganic$Developmental.Stage)


# model
model43_marker_inorganic <- rma.mv(yi=lnRR, 
                                           V=lnRR.sv, 
                                           mods = ~ Biomarker.Category - 1, 
                                           random = list(~1 | References,
                                                         ~1 | Species,
                                                         ~1 | Species_phylo,
                                                         ~1 | Observations), 
                                           R = list(Species_phylo = phylo_cor),
                                           Rscale = "cor",
                                           data = df_inorganic, 
                                           method = "REML", 
                                           #control=list(optimizer="optimParallel",ncpus=2),
                                           sparse = F, 
                                           verbose = T)

#saveRDS(object = model43_marker_inorganic, file = "./results/models/model43_marker_inorganic.RDS")
#model43_marker_inorganic <- readRDS(file = "./results/models/model43_marker_inorganic.RDS")

summary(model43_marker_inorganic)
round(r2_ml(model = model43_marker_inorganic),  # % var explained by moderator
      digits = 4) 


# plot
marker_inorganic_plot <- orchard_plot(object = model43_marker_inorganic, 
                                                  mod = "Biomarker.Category", 
                                          group = "References",
                                          trunk.size = 10,
                                          cb = FALSE,
                                          data = df_inorganic,
                                          xlab = "lnRR", 
                                          transfm = "none") +
  theme(axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) 

ggsave(filename = "./plots/model43_marker_inorganic.pdf", 
       plot = marker_inorganic_plot, 
       device = "pdf", 
       height = 100, 
       width = 200, 
       units = "mm")

##### Model 4.4: biomarkers for Others organic compound #####

# data for model
df_organic <- df00 %>% 
  filter(Pollutant.Class_2b == "Others_organic_compound")
table(df_organic$Developmental.Stage)


# model
model44_marker_organic <- rma.mv(yi=lnRR, 
                                   V=lnRR.sv, 
                                   mods = ~ Biomarker.Category - 1, 
                                   random = list(~1 | References,
                                                 ~1 | Species,
                                                 ~1 | Species_phylo,
                                                 ~1 | Observations), 
                                   R = list(Species_phylo = phylo_cor),
                                   Rscale = "cor",
                                   data = df_organic, 
                                   method = "REML", 
                                   #control=list(optimizer="optimParallel",ncpus=2),
                                   sparse = F, 
                                   verbose = T)

#saveRDS(object = model44_marker_organic, file = "./results/models/model44_marker_organic.RDS")
#model44_marker_organic <- readRDS(file = "./results/models/model44_marker_organic.RDS")

summary(model44_marker_organic)
round(r2_ml(model = model44_marker_organic),  # % var explained by moderator
      digits = 4) 


# plot
marker_organic_plot <- orchard_plot(object = model44_marker_organic, 
                                          mod = "Biomarker.Category", 
                                        group = "References",
                                        trunk.size = 10,
                                        cb = FALSE,
                                        data = df_organic,
                                        xlab = "lnRR", 
                                        transfm = "none") +
  theme(axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) 

ggsave(filename = "./plots/model44_marker_organic.pdf", 
       plot = marker_organic_plot, 
       device = "pdf", 
       height = 100, 
       width = 200, 
       units = "mm")

##### Model 4.5: biomarkers for pesticide #####

# data for model
df_pesticide <- df00 %>% 
  filter(Pollutant.Class_2b == "Pesticide")
table(df_pesticide$Developmental.Stage)


# model
model45_marker_pesticide <- rma.mv(yi=lnRR, 
                                 V=lnRR.sv, 
                                 mods = ~ Biomarker.Category - 1, 
                                 random = list(~1 | References,
                                               ~1 | Species,
                                               ~1 | Species_phylo,
                                               ~1 | Observations), 
                                 R = list(Species_phylo = phylo_cor),
                                 Rscale = "cor",
                                 data = df_pesticide, 
                                 method = "REML", 
                                 #control=list(optimizer="optimParallel",ncpus=2),
                                 sparse = F, 
                                 verbose = T)

#saveRDS(object = model45_marker_pesticide, file = "./results/models/model45_marker_pesticide.RDS")
#model45_marker_pesticide <- readRDS(file = "./results/models/model45_marker_pesticide.RDS")

summary(model45_marker_pesticide)
round(r2_ml(model = model45_marker_pesticide),  # % var explained by moderator
      digits = 4) 


# plot
marker_pesticide_plot <- orchard_plot(object = model45_marker_pesticide, 
                                        mod = "Biomarker.Category", 
                                        group = "References",
                                        trunk.size = 10,
                                        cb = FALSE,
                                        data = df_pesticide,
                                        xlab = "lnRR", 
                                        transfm = "none") +
  theme(axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) 

ggsave(filename = "./plots/model45_marker_pesticide.pdf", 
       plot = marker_pesticide_plot, 
       device = "pdf", 
       height = 75, 
       width = 200, 
       units = "mm")


