##
##
##### Meta-analysis of Colette's dataset
##### Date: Feb 2022
##### PCL
##
##

rm(list = ls())

##
## Script 

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
##### Life mode model #####
##
##

# test datset to test model
df_test <- df00 %>% 
  group_by(Species) %>% 
  filter(lnRR.sv < 10) %>% 
  filter(row_number() == 1) %>% 
  rbind(., 
        df00[sample(x = 1:nrow(df00), 
                    size = 500, 
                    replace = F),])

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
#sparse = T: 
#sparse = F: 
#sparse = F: 
#sparse = F: 

system.time(
  overall_model <- rma.mv(yi=lnRR, 
                          V=lnRR.sv, 
                          mods = ~ Life_cycle2 - 1,
                          random = list(~1 | References,
                                        ~1 | Species,
                                        ~1 | Species_phylo,
                                        ~1 | Observations), 
                          R = list(Species_phylo = phylo_cor),
                          Rscale = "cor",
                          data = df00, 
                          #data = df_test, 
                          method = "REML", 
                          control=list(optimizer="optimParallel",ncpus=3),
                          sparse = F)
)

summary(overall_model) # summary
round(i2_ml(model = overall_model, # heterogeneity values 
            method = "ns"), 
      digits = 4) 

# save model output
#saveRDS(object = overall_model, file = "./results/models/life_cycle_model.RDS")


##
## read models back and plot
overall_model <- readRDS("./results/models/life_cycle_model.RDS")          # warning present

##
## compare outputs
summary(model1)
summary(model1_25)

round(i2_ml(model = model1, method = "ns"), 
      digits = 4) 
round(i2_ml(model = model1_25, method = "ns"), 
      digits = 4) 

##
## EXTREMELLY SIMILAR OUTPUTS, WARNING IS NO PROBLEM IN THIS CASE

##
##
##### Plot Overall model #####
##
##
lifehistory_plot <- orchard_plot(object = overall_model, 
                                 mod = "Life_cycle2", 
                             group = "References",
                             trunk.size = 10,
                             cb = FALSE,
                             data = df00,
                             xlab = "lnRR", 
                             transfm = "none") +
  theme(axis.title = element_text(size = 15),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 15)) +
  scale_fill_manual(values = wes_palette("FantasticFox1", n = 2)) +
  scale_color_manual(values = wes_palette("FantasticFox1", n = 2)) 

ggsave(filename = "./plots/model1_life_history.pdf", 
       plot = lifehistory_plot, 
       device = "pdf", 
       height = 100, 
       width = 200, 
       units = "mm")





