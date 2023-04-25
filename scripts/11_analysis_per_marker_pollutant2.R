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
##### Model 11: Biomarkers and pollutants (2 cat) #####
##
##

##
##### Model 4.1: biomarkers for inorganic #####

# data for model
table(df00$Pollutant.Class_3)
df00$Pollutant.Class_3 <- ifelse(df00$Pollutant.Class_3 == "Organic compound", 
                                 "Organic compound",
                                 "Inorganic compound")

df_inorganic <- df00 %>% 
  filter(Pollutant.Class_3 == "Inorganic compound")


# model
model11_marker_inorganic <- rma.mv(yi=lnRR, 
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
                                   control=list(optimizer="optimParallel",ncpus=3),
                                   sparse = F, 
                                   verbose = T)

saveRDS(object = model11_marker_inorganic, file = "./results/models/model11_marker_inorganic.RDS")
#model11_marker_inorganic <- readRDS(file = "./results/models/model11_marker_inorganic.RDS")

summary(model11_marker_inorganic)
round(r2_ml(model = model11_marker_inorganic),  # % var explained by moderator
      digits = 4) 


# plot
biomarker_inorganic_plot <- orchard_plot(object = model11_marker_inorganic, 
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
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) +
  labs(title = "Inorganic")

ggsave(filename = "./plots/model11_marker_inorganic.pdf", 
       plot = biomarker_inorganic_plot, 
       device = "pdf", 
       height = 125, 
       width = 200, 
       units = "mm")

##
##### Model 4.2: biomarkers for organic #####

# data for model
table(df00$Pollutant.Class_3)
df00$Pollutant.Class_3 <- ifelse(df00$Pollutant.Class_3 == "Organic compound", 
                                 "Organic compound",
                                 "Inorganic compound")

df_organic <- df00 %>% 
  filter(Pollutant.Class_3 == "Organic compound")


# model
model11_marker_organic <- rma.mv(yi=lnRR, 
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
                                   control=list(optimizer="optimParallel",ncpus=3),
                                   sparse = F, 
                                   verbose = T)

saveRDS(object = model11_marker_organic, file = "./results/models/model11_marker_organic.RDS")
#model11_marker_organic <- readRDS(file = "./results/models/model11_marker_organic.RDS")

summary(model11_marker_organic)
round(r2_ml(model = model11_marker_organic),  # % var explained by moderator
      digits = 4) 


# plot
biomarker_organic_plot <- orchard_plot(object = model11_marker_organic, 
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
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) +
  labs(title = "Organic")

ggsave(filename = "./plots/model11_marker_organic.pdf", 
       plot = biomarker_organic_plot, 
       device = "pdf", 
       height = 125, 
       width = 200, 
       units = "mm")
