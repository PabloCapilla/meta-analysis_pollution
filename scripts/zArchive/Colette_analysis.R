##
##
##### Meta-analysis of Colette's dataset
##### Date: Feb 2022
##### PCL
##
##

rm(list = ls())

##
##### libraries ####
##
pacman::p_load(openxlsx, dplyr, tidyr, 
               metafor, ggplot2, extrafont, 
               RColorBrewer, orchaRd, "optimParallel") 
loadfonts()

source("./scripts/0a_R_library/functions.R")
source("./scripts/0a_R_library/orchard_plot_PCL.R")
source("./scripts/0a_R_library/orchard_plot_PCL_noApples.R")

phylo_cor <- readRDS("./results/phylo_cor_full.RDS")

##
##### data #####
##
df00  <- read.xlsx("./data/Final dataset_PB_PCL.xlsx",
                        colNames=T,
                        sheet = 2)
head(df00)

## include observation ID
df00$Observations <- 1:nrow(df00)

##
##
##### Changing species names to match phylo covar matrix and cross-checking #####
##
##
df00[df00$Species == "Lithobates pipiens", "Species"] <- "Rana pipiens" 
df00[df00$Species == "Lithobates catesbeianus", "Species"] <- "Rana catesbeiana" 

#df00$Species <- gsub(x = df00$Species, 
#                             pattern = " ", 
#                             replacement = "_")
df00$Species_phylo <- df00$Species
head(df00)

## matching names between dataset and plylo cor
table(colnames(phylo_cor) %in% df00$Species)
table(colnames(phylo_cor) %in% df00$Species_phylo)

table(df00$Species %in% colnames(phylo_cor))
table(df00$Species_phylo %in% colnames(phylo_cor))


##
##### Calculating adjusted sample sizes #####
##
df00 <- df00 %>% 
  group_by(References, Species, Control.M, Control.SD) %>% 
  mutate(k_shared_control = n()) %>% 
  mutate(Control.N.adj = Control.N / k_shared_control)

## double check this calculation
df00 %>% 
  select(References, 
         Species, 
         Control.M, 
         Control.SD, 
         Control.N,
         k_shared_control,
         Control.N.adj) %>% 
  print(n = 20)

## how many observations n control < 1?
df00 %>%
  ungroup() %>% 
  filter(Control.N.adj < 1) %>% 
  summarise(n_obs = n())

# data for effect size calculation
df01 <- df00 %>%
  filter(Control.N.adj >= 1) 

##
##
##### Effect size calculation #####
##
##

## SMD
df02 <- as.data.frame(escalc(measure = "SMD", 
                             n1i = Control.N.adj, 
                             sd1i = Control.SD, 
                             m1i = Control.M, 
                             n2i = Treatment.N, 
                             sd2i = Treatment.SD, 
                             m2i = Treatment.M, 
                             var.names = c("SMD", "vSMD"), 
                             data = df01,
                             append=TRUE))
head(df02)

## lnRR
df03 <- as.data.frame(escalc(measure = "ROM", 
                             n1i = Control.N.adj, 
                             sd1i = Control.SD, 
                             m1i = Control.M, 
                             n2i = Treatment.N, 
                             sd2i = Treatment.SD, 
                             m2i = Treatment.M, 
                             var.names = c("lnRR","lnRR.sv"), 
                             data = df02,
                             append=TRUE))

## lnCVR
df04 <- as.data.frame(escalc(measure = "CVR", 
                             n1i = Control.N.adj, 
                             sd1i = Control.SD, 
                             m1i = Control.M, 
                             n2i = Treatment.N, 
                             sd2i = Treatment.SD, 
                             m2i = Treatment.M, 
                             var.names = c("lnCVR","lnCVR.sv"), 
                             data = df03,
                             append=TRUE))



##
##
###
##### VALUES WITH SV = O - NEED TO BE CHECKED!!!! #####
###
##
##
df04[df04$lnRR.sv == 0,]

df05 <- df04[df04$lnRR.sv != 0,]

## for publicatoin bias. AMEND ONCE YOU DO THOSE ANALYSES
#data <- data %>% 
#  mutate(sqrt_inv_eff_ss = sqrt((1/Control.N) + (1/Treatment.N)),
#         inv_eff_ss = (1/Control.N) + (1/Treatment.N),
#         publication_year_c = publication_year - mean(publication_year))

##
##
##### Summary of sample sizes #####
##
##
nrow(df05) # total of observations
length(unique(df05$References)) # n of studies
length(unique(df05$Species)) # n of species



##
##### Model 1b: Publicaiton bias #####
##
overall_model_pb <- rma.mv(yi=SMD, 
                           V=vSMD, 
                           mods = ~ 
                             sqrt_inv_eff_ss +
                             publication_year_c ,
                           random = list(~1 | Study.ID,
                                         ~1 | obsID,
                                         ~1 | Species), 
                           data = data, 
                           method = "REML")
summary(overall_model_pb)



##
##### Model 2 : Biomarkers #####
##
biom_model <- rma.mv(yi=SMD, 
                     V=vSMD, 
                     mods = ~ Specific.Biomarker - 1,
                     random = list(~1 | Study.ID,
                                   ~1 | obsID,
                                   ~1 | Species), 
                     data = data, 
                     method = "REML")
summary(biom_model)
round(r2_ml(biom_model), digits = 4) 


biom_plot <- orchard_plot_PCL(object = biom_model, 
                                 mod = "Intercept", 
                                 est_point_size = 4,
                                 alpha = 0.5,
                                 cb = FALSE,
                                 xlab = "SMD", 
                                 ylab = "",
                                 transfm = "none") +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2")

biom_plot <- biom_plot + 
  theme(axis.title = element_text("Arial", size = 10),
        axis.text.x = element_text("Arial", size = 10),
        axis.text.y = element_text("Arial", size = 7, angle = 45),
        legend.position = "none") +
  scale_y_discrete(labels = rev(c("SOD", "MDA", "GSH", "GR", "GPx", "CAT"))) +
  labs(title = "Effect by biomarker")

ggsave(filename = "./plots/biomarker_model.jpeg", 
       plot = biom_plot, 
       device = "jpeg", 
       height = 90, 
       width = 150, 
       units = "mm")


##
##### Model 3 : Pollutant #####
##
pollutant_model <- rma.mv(yi=SMD, 
                     V=vSMD, 
                     mods = ~ Pollutant.Class - 1 ,
                     random = list(~1 | Study.ID,
                                   ~1 | obsID,
                                   ~1 | Species), 
                     data = data, 
                     method = "REML")
summary(pollutant_model)
round(r2_ml(pollutant_model), digits = 4) 


pollutant_plot <- orchard_plot_PCL(object = pollutant_model, 
                              mod = "Intercept", 
                              est_point_size = 4,
                              alpha = 0.5,
                              cb = FALSE,
                              xlab = "SMD", 
                              ylab = "",
                              transfm = "none") +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2")

pollutant_plot <- pollutant_plot + 
  theme(axis.title = element_text("Arial", size = 10),
        axis.text.x = element_text("Arial", size = 10),
        axis.text.y = element_text("Arial", size = 7, angle = 45),
        legend.position = "none") +
  scale_y_discrete(labels = c("Fungicide", "Heavy metal", "Herbicide", 
                              "Insecticide", "Pesticide", "Wastewater contaminant")) +
  labs(title = "Effect by pollutant")

ggsave(filename = "./plots/pollutant_model.jpeg", 
       plot = pollutant_plot, 
       device = "jpeg", 
       height = 90, 
       width = 150, 
       units = "mm")


##
##### Model 4 : Dev stage #####
##
unique(data$Developmental.Stage)
data$Developmental.Stage <- ordered(x = data$Developmental.Stage, 
                                    levels = c("Embryo", "Tadpole", "Adult"))

dev_model <- rma.mv(yi=SMD, 
                    V=vSMD, 
                    mods = ~ Developmental.Stage - 1,
                    random = list(~1 | Study.ID,
                                  ~1 | obsID,
                                  ~1 | Species), 
                    data = data, 
                    method = "REML")
summary(dev_model)
round(r2_ml(dev_model), digits = 4) 


dev_plot <- orchard_plot_PCL(object = dev_model, 
                                   mod = "Intercept", 
                                   est_point_size = 4,
                                   alpha = 0.5,
                                   cb = FALSE,
                                   xlab = "SMD", 
                                   ylab = "",
                                   transfm = "none") +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2")

dev_plot <- dev_plot + 
  theme(axis.title = element_text("Arial", size = 10),
        axis.text.x = element_text("Arial", size = 10),
        axis.text.y = element_text("Arial", size = 7, angle = 45),
        legend.position = "none") +
  scale_y_discrete(labels = c("Embryo", "Tadpole", "Adult")) +
  labs(title = "Effect by developmental stage")

ggsave(filename = "./plots/dev_model.jpeg", 
       plot = dev_plot, 
       device = "jpeg", 
       height = 90, 
       width = 150, 
       units = "mm")


##
##### Model 5 : Levels of Exposure #####
##
unique(data$Level.of.Exposure.T)
data$Level.of.Exposure.T <- ordered(x = data$Level.of.Exposure.T, 
                                    levels = c("Low", "Medium", "High"))

# model below can't be fitted
exp_model <- rma.mv(yi=SMD, 
                     V=vSMD, 
                     mods = ~ 
                       Level.of.Exposure.T - 1 ,
                     random = list(~1 | Study.ID,
                                   ~1 | obsID,
                                   ~1 | Species), 
                     data = data, 
                     method = "REML")
summary(exp_model)
round(r2_ml(exp_model), digits = 4) 


exp_plot <- orchard_plot_PCL(object = exp_model, 
                             mod = "Intercept", 
                             est_point_size = 4,
                             alpha = 0.5,
                             cb = FALSE,
                             xlab = "SMD", 
                             ylab = "",
                             transfm = "none") +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2")

exp_plot <- exp_plot + 
  theme(axis.title = element_text("Arial", size = 10),
        axis.text.x = element_text("Arial", size = 10),
        axis.text.y = element_text("Arial", size = 7, angle = 45),
        legend.position = "none") +
  labs(title = "Effect by Dev Stage")

ggsave(filename = "./plots/exp_model.jpeg", 
       plot = exp_plot, 
       device = "jpeg", 
       height = 90, 
       width = 150, 
       units = "mm")

Specific.Biomarker - 1
Pollutant.Class - 1 
Developmental.Stage - 1 

overall_CVR <- rma.mv(yi=lnCVR, 
                      V=lnCVR.sv, 
                      mods = ~ Pollutant.Class - 1 ,
                      random = list(~1 | Study.ID,
                                    ~1 | obsID,
                                    ~1 | Species), 
                      data = data, 
                      method = "REML")
summary(overall_CVR)
round(i2_ml(overall_CVR), digits = 4) # heterogeneity random effects
round(r2_ml(overall_CVR), digits = 4)

orchard_plot_PCL(object = overall_CVR, 
                 mod = "Intercept", 
                 est_point_size = 4,
                 alpha = 0.5,
                 cb = FALSE,
                 xlab = "SMD", 
                 ylab = "",
                 transfm = "none") +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  theme(axis.title = element_text("Arial", size = 10),
       axis.text.x = element_text("Arial", size = 10),
       axis.text.y = element_text("Arial", size = 7, angle = 45),
       legend.position = "none")


exp_CVR_model <- rma.mv(yi=lnCVR, 
                    V=lnCVR.sv, 
                    mods = ~ 
                      Level.of.Exposure.T - 1 ,
                    random = list(~1 | Study.ID,
                                  ~1 | obsID,
                                  ~1 | Species), 
                    data = data, 
                    method = "REML")
summary(exp_CVR_model)
round(r2_ml(exp_CVR_model), digits = 4) 

orchard_plot_PCL(object = exp_CVR_model, 
                 mod = "Intercept", 
                 est_point_size = 4,
                 alpha = 0.5,
                 cb = FALSE,
                 xlab = "SMD", 
                 ylab = "",
                 transfm = "none") +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2")
