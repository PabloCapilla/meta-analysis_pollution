##
##
##### Meta-analysis of Colette's dataset
##### Date: Feb 2022
##### PCL
##
##

rm(list = ls())

##
## Script 4: meta-analysis per biomaker and dev stage

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
df00$Biomarker.Category_3 <- ifelse(df00$Biomarker.Category == "Indicator", "Indicator", "Antioxidant")
table(df00$Biomarker.Category_3)
table(df00$Biomarker.Category_2)
df00$Biomarker.Category_2 <- df00$Biomarker.Category_3
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
##### Sample size per biomarker and dev stage #####
##
##
table(df00$Developmental.Stage, df00$Biomarker.Category_2)



##
##
##
##### Model 3: Biomarkers and dev stage #####
##
##

##
##### Model 3.1: biomarkers for embryos #####

# data for model
df_embryo <- df00 %>% 
  filter(Developmental.Stage == "Embryo")

# model
model31_marker_embryo <- rma.mv(yi=lnRR, 
                                V=lnRR.sv, 
                                mods = ~ Biomarker.Category_2 - 1,
                                random = list(~1 | References,
                                              ~1 | Species,
                                              ~1 | Species_phylo,
                                              ~1 | Observations), 
                                R = list(Species_phylo = phylo_cor),
                                Rscale = "cor",
                                data = df_embryo, 
                                method = "REML",
                                verbose=TRUE,
                                control=list(optimizer="optimParallel",ncpus=3),
                                sparse = F)
#saveRDS(object = model31_marker_embryo, file = "./results/models/model31_marker_embryo.RDS")
#model31_marker_embryo <- readRDS(file = "./results/models/model31_marker_embryo.RDS")

summary(model31_marker_embryo)
round(r2_ml(model = model31_marker_embryo),  # % var explained by moderator
      digits = 4) 


# plot
biomarker_embryo_plot <- orchard_plot(object = model31_marker_embryo, 
                                      mod = "Biomarker.Category", 
                                      group = "References",
                                      trunk.size = 10,
                                      cb = FALSE,
                                      data = df_embryo,
                                      xlab = "lnRR", 
                                      transfm = "none") +
  theme(axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) 

ggsave(filename = "./plots/model31_biomarker_embryo.pdf", 
       plot = biomarker_embryo_plot, 
       device = "pdf", 
       height = 100, 
       width = 200, 
       units = "mm")

biomarker_embryo_plot2 <- biomarker_embryo_plot +
  scale_x_continuous(limits = c(-3,3))

ggsave(filename = "./plots/model31_biomarker_embryo_STD.pdf", 
       plot = biomarker_embryo_plot2, 
       device = "pdf", 
       height = 75, 
       width = 200, 
       units = "mm")


##
##### Model 3.2: biomarkers for tadpole #####

# data for model
df_tadpole <- df00 %>% 
  filter(Developmental.Stage == "Tadpole")

df_tadpole_trimmed <- df_tadpole %>% 
  filter(lnRR.sv > sort(df_tadpole$lnRR.sv)[10]) %>% 
  filter(lnRR.sv < sort(df_tadpole$lnRR.sv, decreasing = T)[10])

# model
model32_marker_Tadpole <- rma.mv(yi=lnRR, 
                                 V=lnRR.sv, 
                                 mods = ~ Biomarker.Category_2 - 1,
                                 random = list(~1 | References,
                                               ~1 | Species,
                                               ~1 | Species_phylo,
                                               ~1 | Observations), 
                                 R = list(Species_phylo = phylo_cor),
                                 Rscale = "cor",
                                 #data = df_tadpole, 
                                 data = df_tadpole_trimmed, 
                                 method = "REML", 
                                 control=list(optimizer="optimParallel",ncpus=3),
                                 sparse = F)
#saveRDS(object = model32_marker_Tadpole, file = "./results/models/model32_marker_tadpole.RDS")
saveRDS(object = model32_marker_Tadpole, file = "./results/models/model32_marker_tadpole_trimmed_biomarker2.RDS")

##
## read models back and plot
#model32_marker_Tadpole <- readRDS("./results/models/model32_marker_tadpole.RDS")          # warning present
#model32_marker_Tadpole_trimmed <- readRDS("./results/models/model32_marker_tadpole_trimmed.RDS") # warning not present

##
## compare outputs
summary(model32_marker_Tadpole)
summary(model32_marker_Tadpole_trimmed)

round(r2_ml(model = model32_marker_Tadpole), 
      digits = 4) 
round(r2_ml(model = model32_marker_Tadpole_trimmed), 
      digits = 4) 


# plot
biomarker_tadpole_plot <- orchard_plot(object = model32_marker_Tadpole, 
                                           mod = "Biomarker.Category", 
                                           group = "References",
                                           trunk.size = 10,
                                           cb = FALSE,
                                           data = df_tadpole,
                                           xlab = "lnRR", 
                                           transfm = "none") +
  theme(axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) 

ggsave(filename = "./plots/model32_biomarker_tadpole_Biomarker2.pdf", 
       plot = biomarker_tadpole_plot, 
       device = "pdf", 
       height = 75, 
       width = 200, 
       units = "mm")

biomarker_tadpole_plot2 <- biomarker_tadpole_plot +
  scale_x_continuous(limits = c(-3,3))

ggsave(filename = "./plots/model32_biomarker_tadpole_STD.pdf", 
       plot = biomarker_tadpole_plot2, 
       device = "pdf", 
       height = 75, 
       width = 200, 
       units = "mm")
##
##### Model 3.3: biomarkers for adults #####

# data for model
df_adult <- df00 %>% 
  filter(Developmental.Stage == "Adult")

# model
model33_marker_adult <- rma.mv(yi=lnRR, 
                               V=lnRR.sv, 
                               mods = ~ Biomarker.Category_2 - 1,
                               random = list(~1 | References,
                                             ~1 | Species,
                                             ~1 | Species_phylo,
                                             ~1 | Observations), 
                               R = list(Species_phylo = phylo_cor),
                               Rscale = "cor",
                               data = df_adult, 
                               method = "REML", 
                               control=list(optimizer="optimParallel",ncpus=3),
                               sparse = F)
#saveRDS(object = model33_marker_adult, file = "./results/models/model33_marker_adult_Biomarker2.RDS")
#model33_marker_adult <- readRDS("./results/models/model33_marker_adult.RDS")

summary(model33_marker_adult)
round(r2_ml(model = model33_marker_adult),  # % var explained by moderator
      digits = 4) 


# plot
biomarker_adult_plot <- orchard_plot(object = model33_marker_adult, 
                                         mod = "Biomarker.Category", 
                                         group = "References",
                                         trunk.size = 10,
                                         cb = FALSE,
                                         data = df_adult,
                                         xlab = "lnRR", 
                                         transfm = "none") +
  theme(axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) 

ggsave(filename = "./plots/model33_biomarker_adult_biomarker2.pdf", 
       plot = biomarker_adult_plot, 
       device = "pdf", 
       height = 75, 
       width = 200, 
       units = "mm")

biomarker_adult_plot2 <- biomarker_adult_plot +
  scale_x_continuous(limits = c(-3,3))

ggsave(filename = "./plots/model33_biomarker_adult_STD.pdf", 
       plot = biomarker_adult_plot2, 
       device = "pdf", 
       height = 75, 
       width = 200, 
       units = "mm")
