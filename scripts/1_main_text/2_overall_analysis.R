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
#' Overall meta-analysis
#' 
##

##
##### libraries & functions ####
##
pacman::p_load(openxlsx, dplyr, tidyr, metafor, ggplot2, orchaRd) 
source("./scripts/0a_R_library/functions.R")

#####

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

##
## list of papers included in analysis (run only once)
#write.csv(x = df00 %>% 
#            group_by(References) %>% 
#            filter(row_number() == 1) %>% 
#            select(References) %>% 
#            arrange(References), 
#          file = "./list_papers_meta-analysis.csv", 
#          row.names = F)


## matching names between dataset and plylo cor?
table(colnames(phylo_cor) %in% df00$Species) # yes
table(colnames(phylo_cor) %in% df00$Species_phylo) # yes

table(df00$Species %in% colnames(phylo_cor)) # yes
table(df00$Species_phylo %in% colnames(phylo_cor)) # yes

#####

##
##### sample size developmental phase, etc#####
##
table(df00$Developmental.Stage)
table(df00$Developmental.Stage, df00$Biomarker.Category, df00$Pollutant.Class_3) # number of estimates

## number of studies per developmental phase, etc
df00 %>% 
  group_by(Developmental.Stage, Biomarker.Category, Pollutant.Class_3, References) %>% 
  filter(row_number() == 1) %>% 
  group_by(Pollutant.Class_3, Developmental.Stage, Biomarker.Category) %>% 
  summarise(n_studies = n())


#####

##
##### Overall model : MODEL 1 #####
##
##



##
##### Model 1 #####
##
overall_model <- rma.mv(yi=lnRR, 
                        V=lnRR.sv, 
                        mods = ~ 1,
                        random = list(~1 | References,
                                      ~1 | Species,
                                      ~1 | Species_phylo,
                                      ~1 | Observations), 
                        R = list(Species_phylo = phylo_cor),
                        Rscale = "cor",
                        data = df00, 
                        method = "REML", 
                        control=list(optimizer="optimParallel",ncpus=4),
                        sparse = F)

##
## The rma.mv function produces a warning due to ratio of largest to smallest sampling variance extremely large. 
## We have checked that removing the extreme lnRR values that cause the warning and re-running the model produces very similar results.

summary(overall_model) # summary
round(i2_ml(model = overall_model), digits = 4) # heterogeneity values

# save model output
#saveRDS(object = overall_model, file = "./results/models/overall_model1.RDS")

#####

##
## re model back if re-starting the script
#overall_model <- readRDS(file = "./results/models/overall_model1.RDS")


##
##
##### Plot Overall model #####
##
##
overall_plot <- orchard_plot(object = overall_model, 
                             mod = "1", 
                             group = "References",
                             trunk.size = 10,
                             cb = FALSE,
                             xlab = "lnRR", 
                             transfm = "none") +
  theme(axis.title = element_text(size = 15),
        axis.text.x = element_text(size = 15)) +
  scale_fill_manual(values = "#bdbdbd") +
  scale_color_manual(values = "#bdbdbd") 


ggsave(filename = "./plots/Figure 1a.pdf", 
       plot = overall_plot, 
       device = "pdf", 
       height = 100, 
       width = 250, 
       units = "mm")

#####
