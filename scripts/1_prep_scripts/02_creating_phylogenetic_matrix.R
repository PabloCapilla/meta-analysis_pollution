###
###
#' 
#' Script for:
#' TITLE
#' Martin et al. 
#' Preprint: 
#' 
#' Latest update: 2022/11/30
#' 
###
###

# Clear memory to make sure there are not files loaded that could cause problems
rm(list=ls())

##
##
##### Script aim: #####
#' Creating phylogetic correlation matrix to control model for phylogenetic relationship between species
#'

##### 

##
##
##### libraries #####
##
##
pacman::p_load(openxlsx, 
               dplyr, 
               tidyr, 
               metafor, 
               ggplot2, 
               extrafont, 
               RColorBrewer, 
               orchaRd, 
               rotl, 
               ape, 
               curl, 
               treeio) 
loadfonts()
#BiocManager::install("ggtree") # run if ggtree not installed
library(ggtree)

#####

##
##
##### data #####
##
##
df00  <- readRDS("./data/Martin_etal_analysis_dataset.RDS")
head(df00)

#####

##
##
##### creating phylo rel matrix #####
##
##
taxa.corrected <- tnrs_match_names(names = unique(df00$Species))

# check approximate matches: none it seems, good
taxa.corrected[taxa.corrected$approximate_match==TRUE,]

# check synonyms matches:  none it seems, good
taxa.corrected[taxa.corrected$is_synonym==TRUE,]

# check number of matches: none it seems, good
taxa.corrected[taxa.corrected$number_matches>1,]

# retrieving phylogenetic relationships among taxa in the form of a trimmed sub-tree
tree <- tol_induced_subtree(ott_ids = taxa.corrected[["ott_id"]], 
                            label_format = "name")

#####

##
##
##### Tree Plot #####
##
##
plot(tree, type = "fan")
ggtree(tree, layout='circular')

#####

##
##
##### Check tree #####
##
##

# check for the existence of polytomies
is.binary(tree) # there are no polytomies


# checking that the tree includes every species in data table and vv
tree$tip.label <- gsub("_"," ", tree$tip.label)
intersect(as.character(tree$tip.label), as.character(df00$Species))
setdiff(as.character(df00$Species), as.character(tree$tip.label)) #listed in our database but not in the tree
setdiff(as.character(tree$tip.label),as.character(df00$Species)) # listed in the tree but not in our database

##
##
##### computing correlatin matrix and saving #####
##
##

# compute branch lengths of tree
phylo_branch <- compute.brlen(tree, method = "Grafen", power = 1)

# check tree is ultrametric
is.ultrametric(phylo_branch) # TRUE

# matrix to be included in the models
phylo_cor <- vcv(phylo_branch, cor = T)

# 
# save matrix for analyses
saveRDS(phylo_cor, file = "./data/phylogenetic_correlation_matrix.RDS")

#####
