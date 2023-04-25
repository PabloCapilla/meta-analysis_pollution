###
###
#' Meta-analysis of pollutants Colette, Pat & Pablo
#' Script to create phylogenetic relatedness matrix for overall analysis
#' Date: Feb 2022
#' PCL
###
###

##
##### libraries #####
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

#BiocManager::install("ggtree")
#library(ggtree)



# Clear memory to make sure there are not other files around that could cause problems
rm(list=ls())

##
##### data #####
##
df00  <- readRDS("./results/clean_data/clean_analysis_20221011.RDS")
head(df00)


##
##### creating phylo rel matrix #####
##
taxa.corrected <- tnrs_match_names(names = unique(df00$Species))

# check approximate matches: none it seems, good
taxa.corrected[taxa.corrected$approximate_match==TRUE,]

# check synonyms matches: two. Rana pipiens checked; Rana catesbeiana checked
taxa.corrected[taxa.corrected$is_synonym==TRUE,]

# check number of matches: none it seems, good
taxa.corrected[taxa.corrected$number_matches>1,]

# retrieving phylogenetic relationships among taxa in the form of a trimmed sub-tree
tree <- tol_induced_subtree(ott_ids = taxa.corrected[["ott_id"]], 
                            label_format = "name")


##
##### Tree Plot #####
##

plot(tree, type = "fan")
ggtree(tree, layout='circular')



# check for the existence of polytomies
is.binary(tree) # there are no polytomies


# checking that the tree includes every species in data table and vv
tree$tip.label <- gsub("_"," ", tree$tip.label)
intersect(as.character(tree$tip.label), as.character(df00$Species))
setdiff(as.character(df00$Species), as.character(tree$tip.label)) #listed in our database but not in the tree
setdiff(as.character(tree$tip.label),as.character(df00$Species)) # listed in the tree but not in our database


# save the tree
#saveRDS(tree, file = "./data/2_phylo_trees/tree_full.RDS")


# compute branch lengths of tree
phylo_branch <- compute.brlen(tree, method = "Grafen", power = 1)

# check tree is ultrametric
is.ultrametric(phylo_branch) # TRUE

# matrix to be included in the models
phylo_cor <- vcv(phylo_branch, cor = T)
# 
# save matrix for analyses
saveRDS(phylo_cor, file = "./results/clean_data/20221011_data_phylo_cor.RDS")

