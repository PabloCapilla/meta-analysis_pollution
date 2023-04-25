##
##
##### Meta-analysis of Colette's dataset
##### Date: Feb 2022
##### PCL
##
##

rm(list = ls())

##
## Script 3: meta-analysis per biomarker
##
##### libraries & functions ####
##
pacman::p_load(openxlsx, dplyr, tidyr, 
               metafor, ggplot2, extrafont, wesanderson,
               RColorBrewer, orchaRd, optimParallel) 
loadfonts()

source("./scripts/0a_R_library/functions.R")

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
##### Model 2: Biomarkers #####
##
##

# test datset to test model
#df_test <- df00 %>% 
#  filter(lnRR.sv > sort(df00$lnRR.sv)[10]) %>% 
#  filter(lnRR.sv < sort(df00$lnRR.sv, decreasing = T)[10]) %>% 
#  group_by(Species) %>% 
#  filter(row_number() == 1) %>% 
#  rbind(., 
#        df00[sample(x = 1:nrow(df00), 
#                    size = 500, 
#                    replace = F),])

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
##### Model 1 #####
##

# time tests
system.time(
  model2_marker <- rma.mv(yi=lnRR, 
                          V=lnRR.sv, 
                          mods = ~ Biomarker.Category - 1,
                          random = list(~1 | References,
                                        ~1 | Species,
                                        ~1 | Species_phylo,
                                        ~1 | Observations), 
                          R = list(Species_phylo = phylo_cor),
                          Rscale = "cor",
                          #data = df00, 
                          data = df00_trimmed, 
                          method = "REML", 
                          control=list(optimizer="optimParallel",ncpus=3),
                          sparse = F)
)

summary(model2_marker) # summary
round(r2_ml(model = model2_marker),  # % var explained by moderator
      digits = 4) 


# save model output
#saveRDS(object = model2_marker, file = "./results/models/model2_marker.RDS")
#saveRDS(object = model2_marker, file = "./results/models/model2_marker_trimmed.RDS")


##
## read models back and plot
model2_marker <- readRDS("./results/models/model2_marker.RDS")          # warning present
model2_marker_trimmed <- readRDS("./results/models/model2_marker_trimmed.RDS") # warning not present

##
## compare outputs
summary(model2_marker)
summary(model2_marker_trimmed)

round(r2_ml(model = model2_marker), 
      digits = 4) 
round(r2_ml(model = model2_marker_trimmed), 
      digits = 4) 

##
## EXTREMELLY SIMILAR OUTPUTS, WARNING IS NO PROBLEM IN THIS CASE

##
##
##### Plot Overall model #####
##
##
library(wesanderson)
names(wes_palettes)
biomarker_plot <- orchard_plot(object = model2_marker, 
                               data = df00,
                               mod = "Biomarker.Category", 
                               group = "References",
                               trunk.size = 10,
                               cb = FALSE,
                               xlab = "lnRR", 
                               transfm = "none") +
  theme(axis.title = element_text(size = 15),
        axis.text.x = element_text(size = 15)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 3)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3)) 


ggsave(filename = "./plots/model2_biomarker.pdf", 
       plot = biomarker_plot, 
       device = "pdf", 
       height = 100, 
       width = 200, 
       units = "mm")





