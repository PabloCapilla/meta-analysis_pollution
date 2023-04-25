##
##
##### Meta-analysis of Colette's dataset
##### Date: Feb 2022
##### PCL
##
##

rm(list = ls())

##
## Script 7: meta-analysis per specific biomaker and pollutant

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
phylo_cor <- readRDS("./results/clean_data/20220205_data_phylo_cor.RDS")

##
## effect size data
df00  <- readRDS("./results/clean_data/clean_analysis_20220209.RDS")
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

df_enzymatic <- df00 %>% 
  filter(Biomarker.Category == "Enzymatic")

table(df00$Specific.Biomarker, df00$Biomarker.Category)
table(df_enzymatic$Specific.Biomarker, df_enzymatic$Pollutant.Class_2b)


##
##
##
##### Model 7: Specific Biomarkers and pollutants #####
##
##

##
##### Model 7.1: biomarkers for herbicide #####

# data for model
df_herbicide <- df_enzymatic %>% 
  filter(Pollutant.Class_2b == "Herbicide")


# model
model71_marker_herbicide <- rma.mv(yi=lnRR, 
                                   V=lnRR.sv, 
                                   mods = ~ Specific.Biomarker - 1, 
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

#saveRDS(object = model71_marker_herbicide, file = "./results/models/model71_marker_herbicide.RDS")
#model71_marker_herbicide <- readRDS(file = "./results/models/model71_marker_herbicide.RDS")

summary(model71_marker_herbicide)
round(r2_ml(model = model71_marker_herbicide),  # % var explained by moderator
      digits = 4) 


# plot
biomarker_herbicide_plot <- orchard_plot_PCL(object = model71_marker_herbicide, 
                                             mod = "Specific.Biomarker", 
                                             est_point_size = 5,
                                             alpha = 0.5,
                                             cb = FALSE,
                                             xlab = "lnRR", 
                                             ylab = "Biomarker",
                                             transfm = "none", 
                                             angle = 45) +
  theme(axis.title = element_text("Arial", size = 10),
        axis.text.x = element_text("Arial", size = 10)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 5)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 5)) +
  labs(title = "Herbicide")

ggsave(filename = "./plots/model71_marker_herbicide.jpeg", 
       plot = biomarker_herbicide_plot, 
       device = "jpeg", 
       height = 75, 
       width = 200, 
       units = "mm")


##
##### Model 7.2: biomarkers for Metallic elements #####

# data for model
df_Metallic_elements <- df_enzymatic %>% 
  filter(Pollutant.Class_2b == "Metallic elements")


# model
model72_marker_Metallic_elements <- rma.mv(yi=lnRR, 
                                           V=lnRR.sv, 
                                           mods = ~ Specific.Biomarker - 1, 
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

#saveRDS(object = model72_marker_Metallic_elements, file = "./results/models/model72_marker_Metallic_elements.RDS")
#model72_marker_Metallic_elements<- readRDS(file = "./results/models/model72_marker_Metallic_elements.RDS")

summary(model72_marker_Metallic_elements)
round(r2_ml(model = model72_marker_Metallic_elements),  # % var explained by moderator
      digits = 4) 


# plot
marker_Metallic_elements_plot <- orchard_plot_PCL(object = model72_marker_Metallic_elements, 
                                                  mod = "Specific.Biomarker", 
                                                  est_point_size = 5,
                                                  alpha = 0.5,
                                                  cb = FALSE,
                                                  xlab = "lnRR", 
                                                  ylab = "Biomarker",
                                                  transfm = "none", 
                                                  angle = 45) +
  theme(axis.title = element_text("Arial", size = 10),
        axis.text.x = element_text("Arial", size = 10)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 5)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 5)) +
  labs(title = "Metallic elements")

ggsave(filename = "./plots/model72_marker_Metallic_elements.jpeg", 
       plot = marker_Metallic_elements_plot, 
       device = "jpeg", 
       height = 75, 
       width = 200, 
       units = "mm")




##### Model 7.3: biomarkers for Others Inorganic compound #####

# data for model
df_inorganic <- df_enzymatic %>% 
  filter(Pollutant.Class_2b == "Others_Inorganic_compound")


# model
model73_marker_inorganic <- rma.mv(yi=lnRR, 
                                   V=lnRR.sv, 
                                   mods = ~ Specific.Biomarker - 1, 
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

#saveRDS(object = model73_marker_inorganic, file = "./results/models/model73_marker_inorganic.RDS")
#model73_marker_inorganic <- readRDS(file = "./results/models/model73_marker_inorganic.RDS")

summary(model73_marker_inorganic)
round(r2_ml(model = model73_marker_inorganic),  # % var explained by moderator
      digits = 4) 


# plot
marker_inorganic_plot <- orchard_plot_PCL(object = model73_marker_inorganic, 
                                          mod = "Specific.Biomarker", 
                                          est_point_size = 5,
                                          alpha = 0.5,
                                          cb = FALSE,
                                          xlab = "lnRR", 
                                          ylab = "Biomarker",
                                          transfm = "none", 
                                          angle = 45) +
  theme(axis.title = element_text("Arial", size = 10),
        axis.text.x = element_text("Arial", size = 10)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 5)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 5)) +
  labs(title = "Other inorganic")

ggsave(filename = "./plots/model73_marker_inorganic.jpeg", 
       plot = marker_inorganic_plot, 
       device = "jpeg", 
       height = 75, 
       width = 200, 
       units = "mm")

##### Model 7.4: biomarkers for Others organic compound #####

# data for model
df_organic <- df_enzymatic %>% 
  filter(Pollutant.Class_2b == "Others_organic_compound")


# model
model74_marker_organic <- rma.mv(yi=lnRR, 
                                 V=lnRR.sv, 
                                 mods = ~ Specific.Biomarker - 1, 
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

#saveRDS(object = model74_marker_organic, file = "./results/models/model74_marker_organic.RDS")
#model74_marker_organic <- readRDS(file = "./results/models/model74_marker_organic.RDS")

summary(model74_marker_organic)
round(r2_ml(model = model74_marker_organic),  # % var explained by moderator
      digits = 4) 


# plot
marker_organic_plot <- orchard_plot_PCL(object = model74_marker_organic, 
                                        mod = "Specific.Biomarker", 
                                        est_point_size = 5,
                                        alpha = 0.5,
                                        cb = FALSE,
                                        xlab = "lnRR", 
                                        ylab = "Biomarker",
                                        transfm = "none", 
                                        angle = 45) +
  theme(axis.title = element_text("Arial", size = 10),
        axis.text.x = element_text("Arial", size = 10)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 5)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 5)) +
  labs(title = "Other organic")

ggsave(filename = "./plots/model74_marker_organic.jpeg", 
       plot = marker_organic_plot, 
       device = "jpeg", 
       height = 75, 
       width = 200, 
       units = "mm")

##### Model 7.5: biomarkers for pesticide #####

# data for model
df_pesticide <- df_enzymatic %>% 
  filter(Pollutant.Class_2b == "Pesticide")


# model
model75_marker_pesticide <- rma.mv(yi=lnRR, 
                                   V=lnRR.sv, 
                                   mods = ~ Specific.Biomarker - 1, 
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

#saveRDS(object = model75_marker_pesticide, file = "./results/models/model75_marker_pesticide.RDS")
#model75_marker_pesticide <- readRDS(file = "./results/models/model75_marker_pesticide.RDS")

summary(model75_marker_pesticide)
round(r2_ml(model = model75_marker_pesticide),  # % var explained by moderator
      digits = 4) 


# plot
marker_pesticide_plot <- orchard_plot_PCL(object = model75_marker_pesticide, 
                                          mod = "Specific.Biomarker", 
                                          est_point_size = 5,
                                          alpha = 0.5,
                                          cb = FALSE,
                                          xlab = "lnRR", 
                                          ylab = "Biomarker",
                                          transfm = "none", 
                                          angle = 45) +
  theme(axis.title = element_text("Arial", size = 10),
        axis.text.x = element_text("Arial", size = 10)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 5)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 5)) +
  labs(title = "Pesticide")

ggsave(filename = "./plots/model75_marker_pesticide.jpeg", 
       plot = marker_pesticide_plot, 
       device = "jpeg", 
       height = 75, 
       width = 200, 
       units = "mm")


