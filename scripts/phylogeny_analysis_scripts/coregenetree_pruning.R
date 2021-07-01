# Author: Bryden Fields
# Email: bryden.fields@york.ac.uk
# Last updated: 01/07/2021


# Initial set-up ---------

library(ape)
library(phytools)
library(stringr)
library(data.table)

#set the working directory to the directory this script is in.
#setwd("~/what/ever/folder/you/are/working/from") 

# for more info visit:
# http://blog.phytools.org/2011/03/prune-tree-to-list-of-taxa.html
# http://www.phytools.org/Cordoba2017/ex/2/Intro-to-phylogenies.html


# pruning newick tree ------------------

# Read in the 305 core gene newick tree from Cavassim et al., 2020. 
core_tree <- read.tree("../../data/raw_data/phylogeny_analysis_data/305_core_genes.nw")

# view plot of tree
plotTree(core_tree, ftype = "i", fsize = 0.6, lwd = 1)

# how many strains in tree? (should be 196)
Ntip(core_tree)
# [1] 196

# strains in the phylogeny
core_tree$tip.label

# prune with phytools for our 24 strains
# strains we want to keep
ourstrains <- c("-152B", "-147A", "-137B", "-158", "-152A", "-170C", "-145B", "-157B", "-154C", 
                "-165A", "-144A", "-122A", "-149A", "-41", "-126B", "-53", "-135B", "-57", "-135A",
                "-77", "-159", "-74", "-168A", "-67")

# put ourstrain names that match those in the phylogeny into a new vector
pr_strainnames <- core_tree$tip.label %like% ourstrains

pr_strainnames <- core_tree$tip.label[grepl(paste0(ourstrains, collapse = "|"), core_tree$tip.label, ignore.case = TRUE)]
pr_strainnames


# drop all strain tips except our 24 strains
core_tree_subset <- drop.tip(core_tree,
                             setdiff(core_tree$tip.label, pr_strainnames))

# replot tree to check only have 24 strains we want in the tree.
plotTree(core_tree_subset, ftype="i", fsize=0.6, lwd=1)


# save tree (will plot the tree with MEGA)
write.tree(core_tree_subset, file = "../../data/intermediate_data/phylogeny_analysis_data/305_core_genes_24subset.nw")

# make tree in MEGA with newick file.
