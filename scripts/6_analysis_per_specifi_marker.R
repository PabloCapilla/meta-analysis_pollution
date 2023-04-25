##
##
##### Meta-analysis of Colette's dataset
##### Date: Feb 2022
##### PCL
##
##

rm(list = ls())

##
## Script 6: meta-analysis looking at specific biomarker for adults and tadpoles

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
##### Sample size per specific biomarker and biomarker #####
##
##

df_enzymatic <- df00 %>% 
  filter(Biomarker.Category == "Enzymatic")

table(df00$Specific.Biomarker, df00$Biomarker.Category)
table(df_enzymatic$Specific.Biomarker, df_enzymatic$Developmental.Stage)



##
##
##
##### Model 6: Specific biomarkers and dev stage #####
##
##


##
##### Model 6.1: specific biomarkers for tadpoles #####

# data for model
df_tadpole <- df_enzymatic %>% 
  filter(Developmental.Stage == "Tadpole")

# model
model61_spemarker_Tadpole <- rma.mv(yi=lnRR, 
                                    V=lnRR.sv, 
                                    mods = ~ Specific.Biomarker - 1,
                                    random = list(~1 | References,
                                                  ~1 | Species,
                                                  ~1 | Species_phylo,
                                                  ~1 | Observations), 
                                    R = list(Species_phylo = phylo_cor),
                                    Rscale = "cor",
                                    data = df_tadpole, 
                                    #data = df_tadpole_trimmed, 
                                    method = "REML", 
                                    #control=list(optimizer="optimParallel",ncpus=3),
                                    sparse = F,
                                    verbose = T)
#saveRDS(object = model61_spemarker_Tadpole, file = "./results/models/model61_spemarker_Tadpole.RDS")

##
## read models back and plot
#model61_spemarker_Tadpole <- readRDS("./results/models/model32_marker_tadpole.RDS")          # warning present

##
## compare outputs
summary(model61_spemarker_Tadpole)

round(r2_ml(model = model61_spemarker_Tadpole), 
      digits = 4) 

# plot
spebiomarker_tadpole_plot <- orchard_plot(object = model61_spemarker_Tadpole, 
                                           mod = "Specific.Biomarker", 
                                           group = "References",
                                           trunk.size = 10,
                                           cb = FALSE,
                                           data = df_tadpole,
                                           xlab = "lnRR", 
                                           transfm = "none") +
  theme(axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 5)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 5)) 

ggsave(filename = "./plots/model61_spemarker_Tadpole.pdf", 
       plot = spebiomarker_tadpole_plot, 
       device = "pdf", 
       height = 125, 
       width = 200, 
       units = "mm")

##
##### Model 6.2: biomarkers for adults #####

# data for model
df_adult <- df_enzymatic %>% 
  filter(Developmental.Stage == "Adult")

# model
model62_spemarker_adult <- rma.mv(yi=lnRR, 
                                  V=lnRR.sv, 
                                  mods = ~ Specific.Biomarker - 1,
                                  random = list(~1 | References,
                                                ~1 | Species,
                                                ~1 | Species_phylo,
                                                ~1 | Observations), 
                                  R = list(Species_phylo = phylo_cor),
                                  Rscale = "cor",
                                  data = df_adult, 
                                  method = "REML", 
                                  #control=list(optimizer="optimParallel",ncpus=3),
                                  sparse = F)
#saveRDS(object = model62_spemarker_adult, file = "./results/models/model62_spemarker_adult.RDS")
model62_spemarker_adult <- readRDS("./results/models/model62_spemarker_adult.RDS")

summary(model62_spemarker_adult)
round(r2_ml(model = model62_spemarker_adult),  # % var explained by moderator
      digits = 4) 


# plot
spebiomarker_adult_plot <- orchard_plot(object = model62_spemarker_adult, 
                                         mod = "Specific.Biomarker", 
                                         group = "References",
                                         trunk.size = 10,
                                         cb = FALSE,
                                         data = df_adult,
                                         xlab = "lnRR", 
                                         transfm = "none") +
  theme(axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 5)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 5)) 

ggsave(filename = "./plots/model62_spemarker_adult.pdf", 
       plot = spebiomarker_adult_plot, 
       device = "pdf", 
       height = 125, 
       width = 200, 
       units = "mm")

