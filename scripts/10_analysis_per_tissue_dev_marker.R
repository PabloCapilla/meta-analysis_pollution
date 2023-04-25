##
##
##### Meta-analysis of Colette's dataset
##### Date: Feb 2022
##### PCL
##
##

rm(list = ls())

##
## Script 8: meta-analysis per biomarker (combining enz and non-enz)
##
##### libraries & functions ####
##
pacman::p_load(openxlsx, dplyr, tidyr, 
               metafor, ggplot2, extrafont, 
               RColorBrewer, orchaRd, optimParallel) 
loadfonts()

source("./scripts/0a_R_library/functions.R")
source("./scripts/0a_R_library/orchard_plot_PCL.R")
source("./scripts/0a_R_library/orchard_plot_PCL_noApples.R")

##
##### data #####
##

## 
## phylogeny corr matrix
phylo_cor <- readRDS("./results/clean_data/20220225_data_phylo_cor.RDS")

##
## effect size data
df00  <- readRDS("./results/clean_data/clean_analysis_20220225.RDS")
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


df00$Species_phylo[!df00$Species_phylo %in% colnames(phylo_cor)]

##
##
##### Model 10: analysis for main tissues per biomarker #####
##
##

##
## CHECK OBSERVATIONS WITH V LARGE SAMPLING VARIANCE ##
## two observations with very large sampling variance make the model break. 
## I remove then and need to be checked. High SV due to low adj sample size??

df00 %>% 
  filter(lnRR.sv > 10)

## 
## trimmed dataset removing top 10 and bottom 10 obs lnRR.sv observations
df00_trimmed <- df00 %>% 
  filter(lnRR.sv > sort(df00$lnRR.sv)[10]) %>% 
  filter(lnRR.sv < sort(df00$lnRR.sv, decreasing = T)[10])


##
##### Model 10.1 #####
##
table(df00$Biological.Matrix)
df_Whole_tadpole <- df00 %>% 
  filter(Biological.Matrix == "Whole Body") %>% 
  filter(Developmental.Stage == "Tadpole")

# time tests
system.time(
  model10_whole_tad <- rma.mv(yi=lnRR, 
                          V=lnRR.sv, 
                          mods = ~ Biomarker.Category - 1,
                          random = list(~1 | References,
                                        ~1 | Species,
                                        ~1 | Species_phylo,
                                        ~1 | Observations), 
                          R = list(Species_phylo = phylo_cor),
                          Rscale = "cor",
                          data = df_Whole_tadpole, 
                          method = "REML", 
                          control=list(optimizer="optimParallel",ncpus=3),
                          sparse = F)
)

summary(model10_whole_tad) # summary
round(r2_ml(model = model10_whole_tad),  # % var explained by moderator
      digits = 4) 


# save model output
saveRDS(object = model10_whole_tad, file = "./results/models/model10_whole_tad.RDS")


##
## read models back and plot
#model10_whole_tad <- readRDS("./results/models/model10_whole_tad.RDS")          # warning present

##
##
##### Plot Overall model #####
##
##
library(wesanderson)
names(wes_palettes)
whole_tad_plot <- orchard_plot_PCL(object = model10_whole_tad, 
                                   mod = "Biomarker.Category", 
                                   est_point_size = 5,
                                   alpha = 0.5,
                                   cb = FALSE,
                                   xlab = "lnRR", 
                                   ylab = "Biomarker",
                                   transfm = "none", 
                                   angle = 45) +
  theme(axis.title = element_text("Arial", size = 15),
        axis.text.x = element_text("Arial", size = 15)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) +
  labs(title = "Tadpoles")
  

ggsave(filename = "./plots/model10_whole_tad.jpeg", 
       plot = whole_tad_plot, 
       device = "jpeg", 
       height = 100, 
       width = 200, 
       units = "mm")


#
##### Model 10.2 #####
##
df_Whole_embryo <- df00 %>% 
  filter(Biological.Matrix == "Whole Body") %>% 
  filter(Developmental.Stage == "Embryo")

# time tests
system.time(
  model10_whole_embryo <- rma.mv(yi=lnRR, 
                              V=lnRR.sv, 
                              mods = ~ Biomarker.Category - 1,
                              random = list(~1 | References,
                                            ~1 | Species,
                                            ~1 | Species_phylo,
                                            ~1 | Observations), 
                              R = list(Species_phylo = phylo_cor),
                              Rscale = "cor",
                              data = df_Whole_embryo, 
                              method = "REML", 
                              control=list(optimizer="optimParallel",ncpus=3),
                              sparse = F)
)

summary(model10_whole_embryo) # summary
round(r2_ml(model = model10_whole_embryo),  # % var explained by moderator
      digits = 4) 


# save model output
saveRDS(object = model10_whole_embryo, file = "./results/models/model10_whole_embryo.RDS")


##
## read models back and plot
model10_whole_embryo <- readRDS("./results/models/model10_whole_embryo.RDS")          # warning present


library(wesanderson)
names(wes_palettes)
whole_embryo_plot <- orchard_plot_PCL(object = model10_whole_embryo, 
                                   mod = "Biomarker.Category", 
                                   est_point_size = 5,
                                   alpha = 0.5,
                                   cb = FALSE,
                                   xlab = "lnRR", 
                                   ylab = "Biomarker",
                                   transfm = "none", 
                                   angle = 45) +
  theme(axis.title = element_text("Arial", size = 15),
        axis.text.x = element_text("Arial", size = 15)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) +
  labs(title = "Embryo")

ggsave(filename = "./plots/model10_whole_embryo.jpeg", 
       plot = whole_embryo_plot, 
       device = "jpeg", 
       height = 100, 
       width = 200, 
       units = "mm")










##
##### Model 10.3 #####
##
table(df00[df00$Biological.Matrix == "Liver", c("Developmental.Stage")])

table(df00$Biological.Matrix)
df_liver_tadpole <- df00 %>% 
  filter(Biological.Matrix == "Liver") %>% 
  filter(Developmental.Stage == "Tadpole")

# time tests
system.time(
  model10_liver_tad <- rma.mv(yi=lnRR, 
                              V=lnRR.sv, 
                              mods = ~ Biomarker.Category - 1,
                              random = list(~1 | References,
                                            ~1 | Species,
                                            ~1 | Species_phylo,
                                            ~1 | Observations), 
                              R = list(Species_phylo = phylo_cor),
                              Rscale = "cor",
                              data = df_liver_tadpole, 
                              method = "REML", 
                              control=list(optimizer="optimParallel",ncpus=3),
                              sparse = F)
)

summary(model10_liver_tad) # summary
round(r2_ml(model = model10_liver_tad),  # % var explained by moderator
      digits = 4) 


# save model output
saveRDS(object = model10_liver_tad, file = "./results/models/model10_liver_tad.RDS")


##
## read models back and plot
#model10_liver_tad <- readRDS("./results/models/model10_liver_tad.RDS")          # warning present

##
##
##### Plot Overall model #####
##
##
library(wesanderson)
names(wes_palettes)
liver_tad_plot <- orchard_plot_PCL(object = model10_liver_tad, 
                                   mod = "Biomarker.Category", 
                                   est_point_size = 5,
                                   alpha = 0.5,
                                   cb = FALSE,
                                   xlab = "lnRR", 
                                   ylab = "Biomarker",
                                   transfm = "none", 
                                   angle = 45) +
  theme(axis.title = element_text("Arial", size = 15),
        axis.text.x = element_text("Arial", size = 15)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) +
  labs(title = "Tadpoles")


ggsave(filename = "./plots/model10_liver_tad.jpeg", 
       plot = liver_tad_plot, 
       device = "jpeg", 
       height = 100, 
       width = 200, 
       units = "mm")


#
##### Model 10.4 #####
##
df_liver_adults <- df00 %>% 
  filter(Biological.Matrix == "Liver") %>% 
  filter(Developmental.Stage == "Adult")

# time tests
system.time(
  model10_liver_adults <- rma.mv(yi=lnRR, 
                                 V=lnRR.sv, 
                                 mods = ~ Biomarker.Category - 1,
                                 random = list(~1 | References,
                                               ~1 | Species,
                                               ~1 | Species_phylo,
                                               ~1 | Observations), 
                                 R = list(Species_phylo = phylo_cor),
                                 Rscale = "cor",
                                 data = df_liver_adults, 
                                 method = "REML", 
                                 control=list(optimizer="optimParallel",ncpus=3),
                                 sparse = F)
)

summary(model10_liver_adults) # summary
round(r2_ml(model = model10_liver_adults),  # % var explained by moderator
      digits = 4) 


# save model output
saveRDS(object = model10_liver_adults, file = "./results/models/model10_liver_adults.RDS")


##
## read models back and plot
model10_liver_adults <- readRDS("./results/models/model10_liver_adults.RDS")          # warning present


library(wesanderson)
names(wes_palettes)
liver_adults_plot <- orchard_plot_PCL(object = model10_liver_adults, 
                                      mod = "Biomarker.Category", 
                                      est_point_size = 5,
                                      alpha = 0.5,
                                      cb = FALSE,
                                      xlab = "lnRR", 
                                      ylab = "Biomarker",
                                      transfm = "none", 
                                      angle = 45) +
  theme(axis.title = element_text("Arial", size = 15),
        axis.text.x = element_text("Arial", size = 15)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) +
  labs(title = "Adults")

ggsave(filename = "./plots/model10_liver_adults.jpeg", 
       plot = liver_adults_plot, 
       device = "jpeg", 
       height = 100, 
       width = 200, 
       units = "mm")







